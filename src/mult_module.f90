! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module mult_module
! ------------------------------------------------------------------------------
! Code area 2: matrices
! ------------------------------------------------------------------------------

!!****h* Conquest/mult_module
!!  NAME
!!   mult_module
!!  PURPOSE
!!   Collects together various multiplication operations
!!   Also defines parameters relating to the various matrix
!!   multiplications done
!!  USES
!!   basic_types, common, datatypes, GenBlas, GenComms, matrix_comms_module,
!!   matrix_data, matrix_elements_module, matrix_module, maxima_module,
!!   mult_init_module, multiply_module, numbers, primary_module
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   26/04/00
!!  MODIFICATION HISTORY
!!   08/06/2001 dave
!!    Added ROBODoc header and RCS Id and Log tags, and changed
!!    throughout to use GenComms
!!   30/05/2002 dave
!!    Made all long-range (LH or more) matrices allocatable in
!!    main_matrix_multiply and improved headers slightly.  Also added
!!    RCS static object.
!!   12:29, 04/02/2003 drb
!!    Added Trace
!!   13:51, 10/02/2003 drb
!!    Added end_ops and fmmi to deallocate at the end of a run
!!   14:45, 26/02/2003 drb
!!    Added implicit none
!!   10:42, 06/03/2003 drb
!!    Changed prim to bundle in main_matrix_multiply (ifc problem with
!!    aliasing)
!!   14:28, 16/11/2004 dave
!!    Removed TS_TS_T multiply (problematic !)
!!   10:58, 13/02/2006 drb
!!    Complete new style of matrix manipulation which hides matrix
!!    data from rest of code: new routines, variables etc.
!!   2008/02/06 08:25 dave
!!    Changed for output to file not stdout
!!   2011/06/17 ast
!!    nreqs in matrix_transpose is now allocatable
!!  TODO
!!   17/06/2002 dave
!!    Split main_matrix_multiply into sensible small routines (like
!!    getK, getDOmega) Change immi so that the multiplication
!!    initialisation is done by a subroutine and call it many times
!!   2008/05/22 ast
!!    Added timers
!!   2011/03/28 L.Tong
!!    Added matH_dn, matdH_dn, matL_dn, matK_dn, matphi_dn
!!   2011/09/19 L.Tong
!!    Added matM4_dn, matM12_dn
!!   2011/11/24 L.Tong
!!    Added matU_dn, matUT_dn
!!   2012/03/27 L.Tong
!!   - Changed spin implementation
!!   - matL, matK, matH, matdH, matphi, matU, matUT, matM12, matM4,
!!     matLS, matSL are now all arrays of dimension nspin. matX(1)
!!     correspond to spin up channel always and matX(2) corresponds to
!!     spin down channel. This is true even for spin non-polarised
!!     case.
!!   2013/08/23 M.Arita
!!   - Added subroutine symmetrise_matA, which is an extention of 
!!     symmetrise_L and applicable to any sort of matrices
!!   2013/12/02 M.Arita
!!   - Added parameter LS_T_L for XL-BOMD
!!   2014/01/14  lat
!!   - Added parameters matX, matSX, S_X_SX for exx
!!   2014/09/15 18:30 lat
!!    fixed call start/stop_timer to timer_module (not timer_stdlocks_module !)
!!   2016/05/09 10:56 dave
!!    Added helpful output for temp matrix overrun
!!   2016/08/09 18:00 nakata
!!    Added parameters matSatomf, matKEatomf, matNLatomf, 
!!                     matHatomf, matKatomf, matXatomf,
!!                     matSFcoeff, matSFcoeff_tran, matdSFcoeff, matdSFcoeff_e,
!!                     aSa_sCaTr_aSs, sCa_aSs_sSs, aHa_sCaTr_aHs, sCa_aHs_sHs,
!!                     sCaTr_sSs_aSs, aSs_sCa_aSa, sCaTr_sHs_aHs, aHs_sCa_aHa,
!!                     sSs_sSa_sCa, sHs_sHa_sCa,
!!                     aLa_aSa_aLSa, aLa_aHa_aLHa,
!!                     aSs_trans,   aHs_trans,   SFcoeff_trans,
!!                     aSs_pairind, aHs_pairind, SFcoeff_pairind
!!    for atomic-function(atomf)-based calculations.
!!    Changed names from matSC, matCS, SP_trans, SPpairind, SP_PS_H,     H_SP_SP
!!                    to matAP, matPA, AP_trans, APpairind, AP_PA_aHa, aHa_AP_AP
!!    where A = atomic function and P = projecter function for NL.
!!    matU is changed from (sf,nlpf) to (atomf,nlpf) (and its transpose matUT)
!!   2017/03/08 15:00 nakata
!!    Removed PAOP_PS_H, matdS and matdH which are no longer used
!!   2017/12/05 10:21 dave with TM and NW (Mizuho)
!!    Adding multiplications for NA projectors
!!   2018/11/13 17:30 nakata
!!    Changed matS, matT, matTtran, matKE, matNL and matNA to be spin_SF dependent
!!  SOURCE
!!
module mult_module

  use datatypes
  use matrix_module
  use matrix_data,            only: mx_matrices, matrix_pointer
  use global_module,          only: sf, nlpf, paof, atomf, io_lun, nspin, napf
  use GenComms,               only: cq_abort
  use timer_module,           only: start_timer, stop_timer
  use timer_stdclocks_module, only: tmr_std_allocation

  implicit none
  save

  integer(integ), parameter :: aHa_AP_AP = 1   ! type 2
  integer(integ), parameter :: L_S_LS    = 2   ! type 1
  integer(integ), parameter :: L_H_LH    = 3   ! type 1
  integer(integ), parameter :: T_S_TS    = 4   ! type 1
! integer(integ), parameter :: S_TS_S    = 5   ! type 2, maybe do ST_S_S
! integer(integ), parameter :: S_LS_L    = 6   ! type 2, maybe SL_S_L (prune 7?)
  integer(integ), parameter :: SL_S_SLS  = 7   ! type 1
  integer(integ), parameter :: LH_L_SLS  = 8   ! type 1 or 2 depending on H
  integer(integ), parameter :: LS_L_LSL  = 9   ! type 1
  integer(integ), parameter :: HL_S_LSL  = 10  ! type 1 or 2 depending on H
  integer(integ), parameter :: HLS_LS_L  = 11  ! type 2
  integer(integ), parameter :: LSL_SL_L  = 11  ! type 2
  integer(integ), parameter :: SLS_LS_L  = 12  ! type 2
  integer(integ), parameter :: LHL_SL_S  = 13  ! type 2
  integer(integ), parameter :: LSL_SL_H  = 14  ! type 2
  integer(integ), parameter :: AP_PA_aHa = 15  ! type 1 (possibly 2), H(atomf,nlpf ) * H(nlpf ,atomf) = H(atomf,atomf)
  integer(integ), parameter :: TS_T_T    = 16  ! type 2
  integer(integ), parameter :: TH_T_L    = 17  ! type 2 or 1 depending on H
  integer(integ), parameter :: T_H_TH    = 18  ! type 1
  integer(integ), parameter :: T_L_TL    = 19  ! type 1
  integer(integ), parameter :: TL_T_L    = 20  ! type 2
  integer(integ), parameter :: LS_T_L    = 21  ! type 2
  integer(integ), parameter :: S_X_SX    = 22  ! type 1

  ! The indices for atomf-based matrix multiplications
  !     The values  will be set later.
  !     The indices of matrices are demonstrated with lower-case characters: a=atomf, s=sf
  !     C corresponds to SF coefficients.
  !     L corresponds to local subspace regions in the LFD method.
  integer :: aSa_sCaTr_aSs  ! 23, type 1             , S(atomf,atomf)   * C(sf   ,atomf)^T = S(atomf,sf   )
  integer :: sCa_aSs_sSs    ! 24, type 1             , C(sf   ,atomf)   * S(atomf,sf   )   = S(sf   ,sf   )
  integer :: aHa_sCaTr_aHs  ! 25, type 1             , H(atomf,atomf)   * C(sf   ,atomf)^T = H(atomf,sf   )
  integer :: sCa_aHs_sHs    ! 26, type 1             , C(sf   ,atomf)   * H(atomf,sf   )   = H(sf   ,sf   )
  integer :: sCaTr_sSs_aSs  ! 27, type 2             , C(sf   ,atomf)^T * S(sf   ,sf   )   = S(atomf,sf   )
  integer :: aSs_sCa_aSa    ! 28, type 2             , S(atomf,sf   )   * C(sf   ,atomf)   = S(atomf,atomf)
  integer :: sCaTr_sHs_aHs  ! 29, type 2             , C(sf   ,atomf)^T * H(sf   ,sf   )   = H(atomf,sf   )
  integer :: aHs_sCa_aHa    ! 30, type 2             , H(atomf,sf   )   * C(sf   ,atomf)   = H(atomf,atomf)
  integer :: sSs_sSa_sCa    ! 31, type 2             , S(sf   ,sf   )   * S(sf   ,atomf)   = C(sf   ,atomf)
  integer :: sHs_sHa_sCa    ! 32, type 2             , H(sf   ,sf   )   * H(sf   ,atomf)   = C(sf   ,atomf)
  integer :: aLa_aSa_aLSa   ! 33, type 1, for LFD    , L(atomf,atomf)   * S(atomf,atomf)   = Local S(atomf,atomf)
  integer :: aLa_aHa_aLHa   ! 34, type 1, for LFD    , L(atomf,atomf)   * H(atomf,atomf)   = Local H(atomf,atomf)
  integer :: aNA_NAa_aHa ! type 1 for building H: H(atomf, napf) * H(napf, atomf) = H(atomf, atomf)
  integer :: aHa_aNA_aNA ! type 2 for gradients: aHa_NA_NA

  integer(integ), parameter :: mx_mults  = 36

  type(matrix_mult) :: mult(mx_mults)

  integer(integ), parameter :: S_trans   = 1
  integer(integ), parameter :: L_trans   = 2
  integer(integ), parameter :: T_trans   = 3
  integer(integ), parameter :: AP_trans  = 4
  integer(integ), parameter :: LS_trans  = 5
  integer(integ), parameter :: LH_trans  = 6
  integer(integ), parameter :: LSL_trans = 7
  ! The followings will be set later
  integer :: aSs_trans     ! 8
  integer :: aHs_trans     ! 9
  integer :: SFcoeff_trans ! 10
  integer :: aNA_trans     ! 11
  integer :: aNAa_trans    ! 12

  integer(integ), parameter :: mx_trans = 12

  type(pair_data), allocatable, dimension(:,:) :: pairs
  integer, dimension(:), pointer :: Spairind, Lpairind, Tpairind, &
                                    APpairind, LSpairind, LHpairind, &
                                    LSLpairind, &
                                    aSs_pairind, aHs_pairind, SFcoeff_pairind, &
                                    aNApairind, aNAapairind

  type(matrix_trans), dimension(mx_matrices), target :: ltrans
  type(trans_remote) :: gtrans(mx_trans)

  type(matrix_pointer), allocatable, dimension(:) :: mat_p
  integer :: max_matrices, current_matrix
  ! spin independent matrices
  integer, public :: &
       matAP, matPA, mataNA, matNAa
  ! spin_SF dependent matrices
  integer, allocatable, dimension(:), public :: &
       matS, matT, matTtran, matKE, matNL, matNA
  ! spin dependent matrices
  integer, allocatable, dimension(:), public :: &
       matH, matL, matLS, matSL, matK, matphi, matM12, matM4, matU, &
       matUT, matX, matSX
  ! spin independent atomf-based matrices
  integer, public :: &
       matSatomf, matKEatomf, matNLatomf, matNAatomf
  ! spin dependent atomf-based matrices
  integer, allocatable, dimension(:), public :: &
       matHatomf, matKatomf, matXatomf, &
       matSFcoeff, matSFcoeff_tran, matdSFcoeff, matdSFcoeff_e

  ! matrices related to XLBOMD (for DMM, at present)
   integer, allocatable, dimension(:), public :: matXL, matXLvel    ! (nspin)
   integer, allocatable, dimension(:,:), public :: matXL_store      ! (maxiter_Dissipation, nspin)
   integer :: maxiter_Dissipation  

  integer, allocatable, dimension(:), public :: matrix_index, &
                                                trans_index

  ! Parameters for atom distance checking
  real(double) :: InitAtomicDistance_max, InitAtomicDistance_min
  logical :: flag_check_init_atomic_coords
  
  ! Area identification
  integer, parameter, private :: area = 2

!!***

contains

  !!****f* mult_module/immi *
  !!
  !!  NAME
  !!   immi - initialise matrix multiplication
  !!  USAGE
  !!   immi(Partition structure, Primary set structure,
  !!    Covering set structure, processor id)
  !!   immi(parts,prim,gcs,myid)
  !!  PURPOSE
  !!   Initialises all the indexing needed to perform
  !!   matrix multiplication within Conquest
  !!  INPUTS
  !!   type(group_set) :: parts
  !!   type(primary_set) :: prim
  !!   type(cover_set) :: gcs
  !!   integer :: myid
  !!  USES
  !!   basic_types, matrix_data, matrix_elements_module,
  !!   mult_init_module, maxima_module, datatypes
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   26/04/00
  !!  MODIFICATION HISTORY
  !!   08/06/2001 dave
  !!    Added ROBODoc header
  !!   2008/05/22 ast
  !!    Added timers
  !!   2015/07/08 08:01 dave
  !!    Changed the limit on max_matrices to be a user parameter (mainly
  !!    for band output, so that the user who wants to write lots of bands
  !!    can increase the parameter for the output run)
  !!   2016/08/09 21:30 nakata
  !!    Added lines for atomf-based matrices and LFD matrices
  !!   2017/03/08 15:00 nakata
  !!    Removed dSrange, dHrange, PAOPrange and PAOP_PS_H which are no longer used
  !!   2017/12/05 10:24 dave with TM and NW (Mizuho)
  !!    Adding initialisation for NA projector matrices
  !!  SOURCE
  !!
  subroutine immi(parts, prim, gcs, myid, partial)

    use datatypes
    use basic_types
    use matrix_data
    use matrix_elements_module
    use mult_init_module
    use maxima_module,       only: maxpartsproc, maxnabaprocs
    use GenComms,            only: my_barrier, cq_abort
    use matrix_comms_module, only: find_neighbour_procs
    use global_module,       only: area_matrices, mx_temp_matrices, flag_LFD
    use memory_module,       only: reg_alloc_mem, type_int
    use pseudopotential_common, ONLY: flag_neutral_atom_projector

    implicit none

    ! Passed variables
    type(group_set),   target :: parts
    type(primary_set), target :: prim
    type(cover_set),   target :: gcs
    integer                   :: myid
    integer,         optional :: partial

    ! Local variables
    integer      :: stat, i, j
    real(double) :: ra, rc

    ! call start_timer(tmr_std_matrices)
    max_matrices = mx_matrices + mx_temp_matrices ! Leave plenty of room for temporary matrices
    if (.not. allocated(mat_p)) then
       call start_timer(tmr_std_allocation)
       allocate(mat_p(max_matrices), matrix_index(max_matrices), &
                trans_index(max_matrices), STAT=stat)
       if (stat /= 0) &
            call cq_abort("Error allocating matrix indexing ",max_matrices,stat)
       call reg_alloc_mem(area_matrices,2*max_matrices,type_int)
       call stop_timer(tmr_std_allocation)
       do i = 1, max_matrices
          mat_p(i)%length = 0
          nullify(mat_p(i)%matrix)
          matrix_index(i) = 0
          trans_index(i)  = 0
       end do
    end if
    if (.not. allocated(mat)) then
       call start_timer(tmr_std_allocation)
       allocate(mat(maxpartsproc,mx_matrices), STAT=stat)
       if (stat /= 0) &
            call cq_abort("Error allocating matrix type ", &
                          mx_matrices, maxpartsproc)
       call stop_timer(tmr_std_allocation)
    end if
    if(.not. allocated(halo)) then
       call start_timer(tmr_std_allocation)
       allocate(halo(mx_matrices), STAT=stat)
       if (stat /= 0) &
            call cq_abort("Error allocating halo type ", mx_matrices, stat)
       call stop_timer(tmr_std_allocation)
    end if

    ! Set indices of multiplications and trans for atomf-based-matrices
    if (atomf.ne.sf) then
       aSa_sCaTr_aSs = 23
       sCa_aSs_sSs   = 24
       aHa_sCaTr_aHs = 25
       sCa_aHs_sHs   = 26
       sCaTr_sSs_aSs = 27
       aSs_sCa_aSa   = 28
       sCaTr_sHs_aHs = 29
       aHs_sCa_aHa   = 30
       sSs_sSa_sCa   = 31
       sHs_sHa_sCa   = 32
       if (flag_LFD) then
          aLa_aSa_aLSa  = 33
          aLa_aHa_aLHa  = 34
          aNA_NAa_aHa = 35
          aHa_aNA_aNA = 36
       else
          aNA_NAa_aHa = 33
          aHa_aNA_aNA = 34
       endif
       aSs_trans     = 8
       aHs_trans     = 9
       SFcoeff_trans = 10
       aNA_trans = 11
       aNAa_trans = 12
    else
       aNA_NAa_aHa = 23
       aHa_aNA_aNA = 24
       aNA_trans = 8
       aNAa_trans = S_trans
    endif

    mat%sf1_type = sf
    mat%sf2_type = sf
    ! First, initialise matrix indexing for mults and transposes for
    ! the different ranges
    call matrix_ini(parts, prim, gcs, mat(1:prim%groups_on_node,Lrange),    &
                    Lmatind, rcut(Lrange), myid-1, halo(Lrange), ltrans(Lrange))
    call matrix_ini(parts, prim, gcs, mat(1:prim%groups_on_node,Srange),    &
                    Smatind, rcut(Srange), myid-1, halo(Srange), ltrans(Srange))
    if(flag_check_init_atomic_coords) call check_InterAtomicDistances
    call matrix_ini(parts, prim, gcs, mat(1:prim%groups_on_node,Hrange),    &
                    Hmatind, rcut(Hrange), myid-1, halo(Hrange), ltrans(Hrange))
    mat(1:prim%groups_on_node, APrange)%sf1_type = atomf
    mat(1:prim%groups_on_node, APrange)%sf2_type = nlpf
    call matrix_ini(parts, prim, gcs, mat(1:prim%groups_on_node,APrange),   &
                    APmatind, rcut(APrange), myid-1, halo(APrange),         &
                    ltrans(APrange))
    mat(1:prim%groups_on_node,PArange)%sf1_type = nlpf
    mat(1:prim%groups_on_node,PArange)%sf2_type = atomf
    call matrix_ini(parts, prim, gcs, mat(1:prim%groups_on_node,PArange),   &
                    PAmatind, rcut(APrange), myid-1, halo(PArange),         &
                    ltrans(PArange))
    call matrix_ini(parts, prim, gcs, mat(1:prim%groups_on_node,LSrange),   &
                    LSmatind, rcut(LSrange), myid-1, halo(LSrange),         &
                    ltrans(LSrange))
    call matrix_ini(parts, prim, gcs, mat(1:prim%groups_on_node,LHrange),   &
                    LHmatind, rcut(LHrange), myid-1, halo(LHrange),         &
                    ltrans(LHrange))
    call matrix_ini(parts, prim, gcs, mat(1:prim%groups_on_node,HLrange),   &
                    HLmatind, rcut(LHrange), myid-1, halo(HLrange),         &
                    ltrans(HLrange))
    call matrix_ini(parts, prim, gcs, mat(1:prim%groups_on_node,LSLrange),  &
                    LSLmatind, rcut(LSLrange), myid-1, halo(LSLrange),      &
                    ltrans(LSLrange))
    call matrix_ini(parts, prim, gcs, mat(1:prim%groups_on_node,SLSrange),  &
                    SLSmatind, rcut(SLSrange), myid-1, halo(SLSrange),      &
                    ltrans(SLSrange))
    call matrix_ini(parts, prim, gcs, mat(1:prim%groups_on_node,Trange),    &
                    Tmatind, rcut(Trange), myid-1, halo(Trange),            &
                    ltrans(Trange))
    ! Add global transpose for the type 2 mults T?_T_L ?
    call matrix_ini(parts, prim, gcs, mat(1:prim%groups_on_node,TSrange),   &
                    TSmatind, rcut(TSrange), myid-1, halo(TSrange),         &
                    ltrans(TSrange))
    call matrix_ini(parts, prim, gcs, mat(1:prim%groups_on_node,THrange),   &
                    THmatind, rcut(THrange), myid-1, halo(THrange),         &
                    ltrans(THrange))
    call matrix_ini(parts, prim, gcs, mat(1:prim%groups_on_node,TLrange),   &
                    TLmatind, rcut(TLrange), myid-1, halo(TLrange),         &
                    ltrans(TLrange))
    ! Transposes
    call matrix_ini(parts, prim, gcs, mat(1:prim%groups_on_node,LTrrange),  &
                    LTrmatind, rcut(Lrange), myid-1, halo(LTrrange),        &
                    ltrans(LTrrange))
    call matrix_ini(parts, prim, gcs, mat(1:prim%groups_on_node,SLrange),   &
                    SLmatind, rcut(LSrange), myid-1, halo(SLrange),         &
                    ltrans(SLrange))
    call matrix_ini(parts, prim, gcs, mat(1:prim%groups_on_node,TTrrange),  &
                    TTrmatind, rcut(Trange), myid-1, halo(TTrrange),        &
                    ltrans(TTrrange))
    ! Differentials
    call matrix_ini(parts, prim, gcs, mat(1:prim%groups_on_node,SXrange), &
                    SXmatind, rcut(SXrange), myid-1,                      &
                    halo(SXrange), ltrans(SXrange))
    call matrix_ini(parts, prim, gcs, mat(1:prim%groups_on_node,Xrange), &
                    Xmatind, rcut(Xrange), myid-1,                       &
                    halo(Xrange), ltrans(Xrange))
    if (atomf.ne.sf) then
       mat(1:prim%groups_on_node,aSa_range)%sf1_type = atomf
       mat(1:prim%groups_on_node,aSa_range)%sf2_type = atomf
       call matrix_ini(parts, prim, gcs, mat(1:prim%groups_on_node,aSa_range), &
                       aSa_matind, rcut(aSa_range), myid-1,                 &
                       halo(aSa_range), ltrans(aSa_range))
       mat(1:prim%groups_on_node,aHa_range)%sf1_type = atomf
       mat(1:prim%groups_on_node,aHa_range)%sf2_type = atomf
       call matrix_ini(parts, prim, gcs, mat(1:prim%groups_on_node,aHa_range), &
                       aHa_matind, rcut(aHa_range), myid-1,                 &
                       halo(aHa_range), ltrans(aHa_range))
       mat(1:prim%groups_on_node,STr_range)%sf1_type = sf
       mat(1:prim%groups_on_node,STr_range)%sf2_type = sf
       call matrix_ini(parts, prim, gcs, mat(1:prim%groups_on_node,STr_range), &
                       STr_matind, rcut(STr_range), myid-1,                    &
                       halo(STr_range), ltrans(STr_range))
       mat(1:prim%groups_on_node,HTr_range)%sf1_type = sf
       mat(1:prim%groups_on_node,HTr_range)%sf2_type = sf
       call matrix_ini(parts, prim, gcs, mat(1:prim%groups_on_node,HTr_range), &
                       HTr_matind, rcut(HTr_range), myid-1,                    &
                       halo(HTr_range), ltrans(HTr_range))
       mat(1:prim%groups_on_node,aSs_range)%sf1_type = atomf
       mat(1:prim%groups_on_node,aSs_range)%sf2_type = sf
       call matrix_ini(parts, prim, gcs, mat(1:prim%groups_on_node,aSs_range), &
                       aSs_matind, rcut(aSs_range), myid-1,                    &
                       halo(aSs_range), ltrans(aSs_range))
       mat(1:prim%groups_on_node,aHs_range)%sf1_type = atomf
       mat(1:prim%groups_on_node,aHs_range)%sf2_type = sf
       call matrix_ini(parts, prim, gcs, mat(1:prim%groups_on_node,aHs_range), &
                       aHs_matind, rcut(aHs_range), myid-1,                    &
                       halo(aHs_range), ltrans(aHs_range))
       mat(1:prim%groups_on_node,sSa_range)%sf1_type = sf
       mat(1:prim%groups_on_node,sSa_range)%sf2_type = atomf
       call matrix_ini(parts, prim, gcs, mat(1:prim%groups_on_node,sSa_range), &
                       sSa_matind, rcut(sSa_range), myid-1,                    &
                       halo(sSa_range), ltrans(sSa_range))
       mat(1:prim%groups_on_node,sHa_range)%sf1_type = sf
       mat(1:prim%groups_on_node,sHa_range)%sf2_type = atomf
       call matrix_ini(parts, prim, gcs, mat(1:prim%groups_on_node,sHa_range), &
                       sHa_matind, rcut(sHa_range), myid-1,                    &
                       halo(sHa_range), ltrans(sHa_range))
       mat(1:prim%groups_on_node,SFcoeff_range)%sf1_type = sf
       mat(1:prim%groups_on_node,SFcoeff_range)%sf2_type = atomf
       call matrix_ini(parts, prim, gcs, mat(1:prim%groups_on_node,SFcoeff_range), &
                       SFcoeff_matind, rcut(SFcoeff_range), myid-1,                &
                       halo(SFcoeff_range), ltrans(SFcoeff_range))
       mat(1:prim%groups_on_node,SFcoeffTr_range)%sf1_type = atomf
       mat(1:prim%groups_on_node,SFcoeffTr_range)%sf2_type = sf
       call matrix_ini(parts, prim, gcs, mat(1:prim%groups_on_node,SFcoeffTr_range), &
                       SFcoeffTr_matind, rcut(SFcoeffTr_range), myid-1,              &
                       halo(SFcoeffTr_range), ltrans(SFcoeffTr_range))
       if (flag_LFD) then
          mat(1:prim%groups_on_node,LD_range)%sf1_type = atomf
          mat(1:prim%groups_on_node,LD_range)%sf2_type = atomf
          call matrix_ini(parts, prim, gcs, mat(1:prim%groups_on_node,LD_range), &
                          LD_matind, rcut(LD_range), myid-1,                     &
                          halo(LD_range), ltrans(LD_range))
       endif
    endif ! atomf
    ! NA projectors
    if(flag_neutral_atom_projector) then
       mat(1:prim%groups_on_node, aNArange)%sf1_type = atomf
       mat(1:prim%groups_on_node, aNArange)%sf2_type = napf
       call matrix_ini(parts, prim, gcs, mat(1:prim%groups_on_node,aNArange),   &
            aNAmatind, rcut(aNArange), myid-1, halo(aNArange),         &
            ltrans(aNArange))
       mat(1:prim%groups_on_node,NAarange)%sf1_type = napf
       mat(1:prim%groups_on_node,NAarange)%sf2_type = atomf
       call matrix_ini(parts, prim, gcs, mat(1:prim%groups_on_node,NAarange),   &
            NAamatind, rcut(NAarange), myid-1, halo(NAarange),         &
            ltrans(NAarange))
    end if
    call associate_matrices
    call find_neighbour_procs(parts, halo(max_range))
    call start_timer(tmr_std_allocation)
    allocate(pairs(maxnabaprocs+1,mx_trans), STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error allocating pairs in immi: ", &
                       maxnabaprocs+1, mx_trans)
    call stop_timer(tmr_std_allocation)
    do i = 1, mx_trans
       do j = 1, maxnabaprocs+1
          nullify(pairs(j,i)%submat)
       end do
    end do
    call trans_ini(parts, prim, gcs, mat(1:prim%groups_on_node,Lrange),      &
                   myid-1, halo(Lrange), halo(Lrange), ltrans(Lrange),       &
                   gtrans(L_trans), pairs(:,L_trans), Lpairind)
    call trans_ini(parts, prim, gcs, mat(1:prim%groups_on_node,Srange),      &
                   myid-1, halo(Srange), halo(Srange), ltrans(Srange),       &
                   gtrans(S_trans), pairs(:,S_trans), Spairind)
    call trans_ini(parts, prim, gcs, mat(1:prim%groups_on_node,APrange),     &
                   myid-1, halo(APrange), halo(PArange), ltrans(APrange),    &
                   gtrans(AP_trans), pairs(:,AP_trans), APpairind)
    call trans_ini(parts, prim, gcs, mat(1:prim%groups_on_node,LSrange),     &
                   myid-1, halo(LSrange), halo(LSrange), ltrans(LSrange),    &
                   gtrans(LS_trans), pairs(:,LS_trans), LSpairind)
    call trans_ini(parts, prim, gcs, mat(1:prim%groups_on_node,LHrange),     &
                   myid-1, halo(LHrange), halo(LHrange), ltrans(LHrange),    &
                   gtrans(LH_trans), pairs(:,LH_trans), LHpairind)
    call trans_ini(parts, prim, gcs, mat(1:prim%groups_on_node,LSLrange),    &
                   myid-1, halo(LSLrange), halo(LSLrange), ltrans(LSLrange), &
                   gtrans(LSL_trans), pairs(:, LSL_trans), LSLpairind)
    call trans_ini(parts, prim, gcs, mat(1:prim%groups_on_node,Trange),      &
                   myid-1, halo(Trange), halo(Trange), ltrans(Trange),       &
                   gtrans(T_trans), pairs(:, T_trans), Tpairind)
    if (atomf.ne.sf) then
       call trans_ini(parts, prim, gcs, mat(1:prim%groups_on_node,aSs_range),                &
                      myid-1, halo(aSs_range), halo(sSa_range), ltrans(aSs_range), &
                      gtrans(aSs_trans), pairs(:, aSs_trans), aSs_pairind)
       call trans_ini(parts, prim, gcs, mat(1:prim%groups_on_node,aHs_range),                &
                      myid-1, halo(aHs_range), halo(sHa_range), ltrans(aHs_range), &
                      gtrans(aHs_trans), pairs(:, aHs_trans), aHs_pairind)
       call trans_ini(parts, prim, gcs, mat(1:prim%groups_on_node,SFcoeff_range),                &
                      myid-1, halo(SFcoeff_range), halo(SFcoeffTr_range), ltrans(SFcoeff_range), &
                      gtrans(SFcoeff_trans), pairs(:, SFcoeff_trans), SFcoeff_pairind)
       if(flag_neutral_atom_projector) then
          call trans_ini(parts, prim, gcs, mat(1:prim%groups_on_node,aNArange),     &
               myid-1, halo(aNArange), halo(NAarange), ltrans(aNArange),    &
               gtrans(aNA_trans), pairs(:,aNA_trans), aNApairind)
          call trans_ini(parts, prim, gcs, mat(1:prim%groups_on_node,aHa_range),     &
               myid-1, halo(aHa_range), halo(aHa_range), ltrans(aHa_range),    &
               gtrans(aNAa_trans), pairs(:,aNAa_trans), aNAapairind)
       end if
    ! NA projectors
    else if(flag_neutral_atom_projector) then
       call trans_ini(parts, prim, gcs, mat(1:prim%groups_on_node,aNArange),     &
            myid-1, halo(aNArange), halo(NAarange), ltrans(aNArange),    &
            gtrans(aNA_trans), pairs(:,aNA_trans), aNApairind)
       !aNAapairind = Spairind
       aNAa_trans = S_trans
    end if

    ! Now initialise the matrix multiplications
    ! Somewhere, we need to set mult_type 1/2 - probably here !
    ! 1-5
    ! (a) We need to arrange a global transpose before a type 2 mult
    ! (b) mult_wrap is passed A,B,C  - and send mat_mult C,B,A for a
    !     type 2 mult
    ! (c) HotInvS_mm does TS_T_S as type two, and sends mult_wrap
    !     TS, T, S
    mult(aHa_AP_AP)%mult_type = 2
    mult(aHa_AP_AP)%amat    => mat(1:prim%groups_on_node,APrange)
    mult(aHa_AP_AP)%bmat    => mat(1:prim%groups_on_node,PArange)
    mult(aHa_AP_AP)%cmat    => mat(1:prim%groups_on_node,aHa_range)
    mult(aHa_AP_AP)%ahalo   => halo(APrange)
    mult(aHa_AP_AP)%chalo   => halo(aHa_range)
    mult(aHa_AP_AP)%ltrans  => ltrans(APrange)
    !mult(aHa_AP_AP)%bindex => PAmatind
    mult(aHa_AP_AP)%parts   => parts
    mult(aHa_AP_AP)%prim    => prim
    mult(aHa_AP_AP)%gcs     => gcs
    call mult_ini(mult(aHa_AP_AP), PAmatind, myid-1, prim%n_prim, parts)
    mult(L_S_LS)%mult_type = 1
    mult(L_S_LS)%amat    => mat(1:prim%groups_on_node,Lrange)
    mult(L_S_LS)%bmat    => mat(1:prim%groups_on_node,Srange)
    mult(L_S_LS)%cmat    => mat(1:prim%groups_on_node,LSrange)
    mult(L_S_LS)%ahalo   => halo(Lrange)
    mult(L_S_LS)%chalo   => halo(LSrange)
    mult(L_S_LS)%ltrans  => ltrans(Lrange)
    !mult(L_S_LS)%bindex => Smatind
    mult(L_S_LS)%parts   => parts
    mult(L_S_LS)%prim    => prim
    mult(L_S_LS)%gcs     => gcs
    call mult_ini(mult(L_S_LS), Smatind, myid-1, prim%n_prim, parts)
    mult(L_H_LH)%mult_type = 1
    mult(L_H_LH)%amat    => mat(1:prim%groups_on_node,Lrange)
    mult(L_H_LH)%bmat    => mat(1:prim%groups_on_node,Hrange)
    mult(L_H_LH)%cmat    => mat(1:prim%groups_on_node,LHrange)
    mult(L_H_LH)%ahalo   => halo(Lrange)
    mult(L_H_LH)%chalo   => halo(LHrange)
    mult(L_H_LH)%ltrans  => ltrans(Lrange)
    !mult(L_H_LH)%bindex => Hmatind
    mult(L_H_LH)%parts   => parts
    mult(L_H_LH)%prim    => prim
    mult(L_H_LH)%gcs     => gcs
    call mult_ini(mult(L_H_LH), Hmatind, myid-1, prim%n_prim, parts)
    mult(T_S_TS)%mult_type = 1
    mult(T_S_TS)%amat    => mat(1:prim%groups_on_node,Trange)
    mult(T_S_TS)%bmat    => mat(1:prim%groups_on_node,Srange)
    mult(T_S_TS)%cmat    => mat(1:prim%groups_on_node,TSrange)
    mult(T_S_TS)%ahalo   => halo(Trange)
    mult(T_S_TS)%chalo   => halo(TSrange)
    mult(T_S_TS)%ltrans  => ltrans(Trange)
    !mult(T_S_TS)%bindex => Smatind
    mult(T_S_TS)%parts   => parts
    mult(T_S_TS)%prim    => prim
    mult(T_S_TS)%gcs     => gcs
    call mult_ini(mult(T_S_TS), Smatind, myid-1, prim%n_prim, parts)
    ! Removed 2006/03/28 DRB
    ! Not used: if re-instated, will require ST index
    !mult(S_TS_S)%mult_type = 2
    !mult(S_TS_S)%amat   => mat(1:prim%groups_on_node,Srange)
    !mult(S_TS_S)%bmat   => mat(1:prim%groups_on_node,TSrange)
    !mult(S_TS_S)%cmat   => mat(1:prim%groups_on_node,Srange)
    !mult(S_TS_S)%ahalo  => halo(Srange)
    !mult(S_TS_S)%chalo  => halo(Srange)
    !mult(S_TS_S)%ltrans => ltrans(Srange)
    !mult(S_TS_S)%bindex => TSmatind
    !mult(S_TS_S)%parts  => parts
    !mult(S_TS_S)%prim   => prim
    !mult(S_TS_S)%gcs    => gcs
    !call mult_ini(mult(S_TS_S),myid-1,prim%n_prim,parts)
    ! 6-10
    !mult(S_LS_L)%mult_type = 2
    !mult(S_LS_L)%amat   => mat(1:prim%groups_on_node,Lrange)
    !mult(S_LS_L)%bmat   => mat(1:prim%groups_on_node,LSrange)
    !mult(S_LS_L)%cmat   => mat(1:prim%groups_on_node,Srange)
    !mult(S_LS_L)%ahalo  => halo(Lrange)
    !mult(S_LS_L)%chalo  => halo(Srange)
    !mult(S_LS_L)%ltrans => ltrans(Lrange)
    !mult(S_LS_L)%bindex => LSmatind
    !mult(S_LS_L)%parts  => parts
    !mult(S_LS_L)%prim   => prim
    !mult(S_LS_L)%gcs    => gcs
    !call mult_ini(mult(S_LS_L),myid-1,prim%n_prim,parts)
    mult(SL_S_SLS)%mult_type = 1
    mult(SL_S_SLS)%amat    => mat(1:prim%groups_on_node,SLrange)
    mult(SL_S_SLS)%bmat    => mat(1:prim%groups_on_node,Srange)
    mult(SL_S_SLS)%cmat    => mat(1:prim%groups_on_node,SLSrange)
    mult(SL_S_SLS)%ahalo   => halo(SLrange)
    mult(SL_S_SLS)%chalo   => halo(SLSrange)
    mult(SL_S_SLS)%ltrans  => ltrans(SLrange)
    !mult(SL_S_SLS)%bindex => Smatind
    mult(SL_S_SLS)%parts   => parts
    mult(SL_S_SLS)%prim    => prim
    mult(SL_S_SLS)%gcs     => gcs
    call mult_ini(mult(SL_S_SLS), Smatind, myid-1, prim%n_prim, parts)
    rc = (2.0_double * rcut(Srange) + rcut(Lrange))
    ra = rcut(Lrange) + rcut(Hrange)
    if (rc >= ra) then
       mult(LH_L_SLS)%mult_type = 1
       mult(LH_L_SLS)%amat    => mat(1:prim%groups_on_node,LHrange)
       mult(LH_L_SLS)%bmat    => mat(1:prim%groups_on_node,Lrange)
       mult(LH_L_SLS)%cmat    => mat(1:prim%groups_on_node,SLSrange)
       mult(LH_L_SLS)%ahalo   => halo(LHrange)
       mult(LH_L_SLS)%chalo   => halo(SLSrange)
       mult(LH_L_SLS)%ltrans  => ltrans(LHrange)
       !mult(LH_L_SLS)%bindex => Lmatind
       mult(LH_L_SLS)%parts   => parts
       mult(LH_L_SLS)%prim    => prim
       mult(LH_L_SLS)%gcs     => gcs
       call mult_ini(mult(LH_L_SLS), Lmatind, myid-1, prim%n_prim, parts)
    else
       mult(LH_L_SLS)%mult_type = 2
       mult(LH_L_SLS)%amat    => mat(1:prim%groups_on_node,SLSrange)
       mult(LH_L_SLS)%bmat    => mat(1:prim%groups_on_node,LTrrange)
       mult(LH_L_SLS)%cmat    => mat(1:prim%groups_on_node,LHrange)
       mult(LH_L_SLS)%ahalo   => halo(SLSrange)
       mult(LH_L_SLS)%chalo   => halo(LHrange)
       mult(LH_L_SLS)%ltrans  => ltrans(SLSrange)
       !mult(LH_L_SLS)%bindex => LTrmatind
       mult(LH_L_SLS)%parts   => parts
       mult(LH_L_SLS)%prim    => prim
       mult(LH_L_SLS)%gcs     => gcs
       call mult_ini(mult(LH_L_SLS), LTrmatind, myid-1, prim%n_prim, parts)
    endif
    mult(LS_L_LSL)%mult_type = 1
    mult(LS_L_LSL)%amat    => mat(1:prim%groups_on_node,LSrange)
    mult(LS_L_LSL)%bmat    => mat(1:prim%groups_on_node,Lrange)
    mult(LS_L_LSL)%cmat    => mat(1:prim%groups_on_node,LSLrange)
    mult(LS_L_LSL)%ahalo   => halo(LSrange)
    mult(LS_L_LSL)%chalo   => halo(LSLrange)
    mult(LS_L_LSL)%ltrans  => ltrans(LSrange)
    !mult(LS_L_LSL)%bindex => Lmatind
    mult(LS_L_LSL)%parts   => parts
    mult(LS_L_LSL)%prim    => prim
    mult(LS_L_LSL)%gcs     => gcs
    call mult_ini(mult(LS_L_LSL), Lmatind, myid-1, prim%n_prim, parts)
    !ra = (rcut(Hrange)+rcut(Lrange))
    !rc = 2.0_double*rcut(Lrange)+rcut(Srange)
    !if(rc>=ra) then
       mult(HL_S_LSL)%mult_type = 1
       mult(HL_S_LSL)%amat    => mat(1:prim%groups_on_node,HLrange)
       mult(HL_S_LSL)%bmat    => mat(1:prim%groups_on_node,Srange)
       mult(HL_S_LSL)%cmat    => mat(1:prim%groups_on_node,LSLrange)
       mult(HL_S_LSL)%ahalo   => halo(HLrange)
       mult(HL_S_LSL)%chalo   => halo(LSLrange)
       mult(HL_S_LSL)%ltrans  => ltrans(HLrange)
       !mult(HL_S_LSL)%bindex => Smatind
       mult(HL_S_LSL)%parts   => parts
       mult(HL_S_LSL)%prim    => prim
       mult(HL_S_LSL)%gcs     => gcs
       call mult_ini(mult(HL_S_LSL), Smatind, myid-1, prim%n_prim, parts)
    !else
    !   mult(HL_S_LSL)%mult_type = 2
    !   mult(HL_S_LSL)%amat   => mat(1:prim%groups_on_node,LSLrange)
    !   mult(HL_S_LSL)%bmat   => mat(1:prim%groups_on_node,Srange)
    !   mult(HL_S_LSL)%cmat   => mat(1:prim%groups_on_node,LHrange)
    !   mult(HL_S_LSL)%ahalo  => halo(LSLrange)
    !   mult(HL_S_LSL)%chalo  => halo(LHrange)
    !   mult(HL_S_LSL)%ltrans => ltrans(LSLrange)
    !   mult(HL_S_LSL)%bindex => Smatind
    !   mult(HL_S_LSL)%parts  => parts
    !   mult(HL_S_LSL)%prim   => prim
    !   mult(HL_S_LSL)%gcs    => gcs
    !   call mult_ini(mult(HL_S_LSL), myid-1, prim%n_prim, parts)
    !end if
    ! 11-15
    mult(HLS_LS_L)%mult_type = 2
    mult(HLS_LS_L)%amat    => mat(1:prim%groups_on_node,Lrange)
    mult(HLS_LS_L)%bmat    => mat(1:prim%groups_on_node,SLrange)
    mult(HLS_LS_L)%cmat    => mat(1:prim%groups_on_node,LSLrange)
    mult(HLS_LS_L)%ahalo   => halo(Lrange)
    mult(HLS_LS_L)%chalo   => halo(LSLrange)
    mult(HLS_LS_L)%ltrans  => ltrans(Lrange)
    !mult(HLS_LS_L)%bindex => SLmatind
    mult(HLS_LS_L)%parts   => parts
    mult(HLS_LS_L)%prim    => prim
    mult(HLS_LS_L)%gcs     => gcs
    call mult_ini(mult(HLS_LS_L), SLmatind, myid-1, prim%n_prim, parts)
    mult(SLS_LS_L)%mult_type = 2
    mult(SLS_LS_L)%amat    => mat(1:prim%groups_on_node,Lrange)
    mult(SLS_LS_L)%bmat    => mat(1:prim%groups_on_node,SLrange)
    mult(SLS_LS_L)%cmat    => mat(1:prim%groups_on_node,SLSrange)
    mult(SLS_LS_L)%ahalo   => halo(Lrange)
    mult(SLS_LS_L)%chalo   => halo(SLSrange)
    mult(SLS_LS_L)%ltrans  => ltrans(Lrange)
    !mult(SLS_LS_L)%bindex => SLmatind
    mult(SLS_LS_L)%parts   => parts
    mult(SLS_LS_L)%prim    => prim
    mult(SLS_LS_L)%gcs     => gcs
    call mult_ini(mult(SLS_LS_L), SLmatind, myid-1, prim%n_prim, parts)
    mult(LHL_SL_S)%mult_type = 2
    mult(LHL_SL_S)%amat    => mat(1:prim%groups_on_node,Srange)
    mult(LHL_SL_S)%bmat    => mat(1:prim%groups_on_node,LSrange)
    mult(LHL_SL_S)%cmat    => mat(1:prim%groups_on_node,SLSrange)
    mult(LHL_SL_S)%ahalo   => halo(Srange)
    mult(LHL_SL_S)%chalo   => halo(SLSrange)
    mult(LHL_SL_S)%ltrans  => ltrans(Srange)
    !mult(LHL_SL_S)%bindex => LSmatind
    mult(LHL_SL_S)%parts   => parts
    mult(LHL_SL_S)%prim    => prim
    mult(LHL_SL_S)%gcs     => gcs
    call mult_ini(mult(LHL_SL_S), LSmatind, myid-1, prim%n_prim, parts)
    mult(LSL_SL_H)%mult_type = 2
    mult(LSL_SL_H)%amat    => mat(1:prim%groups_on_node,Hrange)
    mult(LSL_SL_H)%bmat    => mat(1:prim%groups_on_node,LSrange)
    mult(LSL_SL_H)%cmat    => mat(1:prim%groups_on_node,LSLrange)
    mult(LSL_SL_H)%ahalo   => halo(Hrange)
    mult(LSL_SL_H)%chalo   => halo(LSLrange)
    mult(LSL_SL_H)%ltrans  => ltrans(Hrange)
    !mult(LSL_SL_H)%bindex => LSmatind
    mult(LSL_SL_H)%parts   => parts
    mult(LSL_SL_H)%prim    => prim
    mult(LSL_SL_H)%gcs     => gcs
    call mult_ini(mult(LSL_SL_H), LSmatind, myid-1, prim%n_prim, parts)
    ra = rcut(APrange)
    rc = rcut(aHa_range)
    if (rc >= ra) then
       mult(AP_PA_aHa)%mult_type = 1
       mult(AP_PA_aHa)%amat    => mat(1:prim%groups_on_node,APrange)
       mult(AP_PA_aHa)%bmat    => mat(1:prim%groups_on_node,PArange)
       mult(AP_PA_aHa)%cmat    => mat(1:prim%groups_on_node,aHa_range)
       mult(AP_PA_aHa)%ahalo   => halo(APrange)
       mult(AP_PA_aHa)%chalo   => halo(aHa_range)
       mult(AP_PA_aHa)%ltrans  => ltrans(APrange)
       !mult(AP_PA_aHa)%bindex => PAmatind
       mult(AP_PA_aHa)%parts   => parts
       mult(AP_PA_aHa)%prim    => prim
       mult(AP_PA_aHa)%gcs     => gcs
       call mult_ini(mult(AP_PA_aHa), PAmatind, myid-1, prim%n_prim, parts)
    else
       mult(AP_PA_aHa)%mult_type = 2
       mult(AP_PA_aHa)%amat    => mat(1:prim%groups_on_node,aHa_range)
       mult(AP_PA_aHa)%bmat    => mat(1:prim%groups_on_node,APrange)
       mult(AP_PA_aHa)%cmat    => mat(1:prim%groups_on_node,APrange)
       mult(AP_PA_aHa)%ahalo   => halo(aHa_range)
       mult(AP_PA_aHa)%chalo   => halo(APrange)
       mult(AP_PA_aHa)%ltrans  => ltrans(aHa_range)
       !mult(AP_PA_aHa)%bindex => APmatind
       mult(AP_PA_aHa)%parts   => parts
       mult(AP_PA_aHa)%prim    => prim
       mult(AP_PA_aHa)%gcs     => gcs
       call mult_ini(mult(AP_PA_aHa), APmatind, myid-1, prim%n_prim, parts)
    end if
    ! 16-20
    mult(TS_T_T)%mult_type = 2
    mult(TS_T_T)%amat    => mat(1:prim%groups_on_node,Trange)
    mult(TS_T_T)%bmat    => mat(1:prim%groups_on_node,TTrrange)
    mult(TS_T_T)%cmat    => mat(1:prim%groups_on_node,TSrange)
    mult(TS_T_T)%ahalo   => halo(Trange)
    mult(TS_T_T)%chalo   => halo(TSrange)
    mult(TS_T_T)%ltrans  => ltrans(TTrrange)
    !mult(TS_T_T)%bindex => TTrmatind
    mult(TS_T_T)%parts   => parts
    mult(TS_T_T)%prim    => prim
    mult(TS_T_T)%gcs     => gcs
    call mult_ini(mult(TS_T_T), TTrmatind, myid-1, prim%n_prim, parts)
    ra = rcut(THrange)
    rc = rcut(Lrange)
    if (rc >= ra) then
       mult(TH_T_L)%mult_type = 1
       mult(TH_T_L)%amat    => mat(1:prim%groups_on_node,THrange)
       mult(TH_T_L)%bmat    => mat(1:prim%groups_on_node,Trange)
       mult(TH_T_L)%cmat    => mat(1:prim%groups_on_node,Lrange)
       mult(TH_T_L)%ahalo   => halo(THrange)
       mult(TH_T_L)%chalo   => halo(Lrange)
       mult(TH_T_L)%ltrans  => ltrans(THrange)
       !mult(TH_T_L)%bindex => Tmatind
       mult(TH_T_L)%parts   => parts
       mult(TH_T_L)%prim    => prim
       mult(TH_T_L)%gcs     => gcs
       call mult_ini(mult(TH_T_L), Tmatind, myid-1, prim%n_prim, parts)
    else
       mult(TH_T_L)%mult_type = 2
       mult(TH_T_L)%amat    => mat(1:prim%groups_on_node,Lrange)
       mult(TH_T_L)%bmat    => mat(1:prim%groups_on_node,TTrrange)
       mult(TH_T_L)%cmat    => mat(1:prim%groups_on_node,THrange)
       mult(TH_T_L)%ahalo   => halo(Lrange)
       mult(TH_T_L)%chalo   => halo(THrange)
       mult(TH_T_L)%ltrans  => ltrans(Lrange)
       !mult(TH_T_L)%bindex => TTrmatind
       mult(TH_T_L)%parts   => parts
       mult(TH_T_L)%prim    => prim
       mult(TH_T_L)%gcs     => gcs
       call mult_ini(mult(TH_T_L), TTrmatind, myid-1, prim%n_prim, parts)
    end if
    mult(T_H_TH)%mult_type = 1
    mult(T_H_TH)%amat    => mat(1:prim%groups_on_node,Trange)
    mult(T_H_TH)%bmat    => mat(1:prim%groups_on_node,Hrange)
    mult(T_H_TH)%cmat    => mat(1:prim%groups_on_node,THrange)
    mult(T_H_TH)%ahalo   => halo(Trange)
    mult(T_H_TH)%chalo   => halo(THrange)
    mult(T_H_TH)%ltrans  => ltrans(Trange)
    !mult(T_H_TH)%bindex => Hmatind
    mult(T_H_TH)%parts   => parts
    mult(T_H_TH)%prim    => prim
    mult(T_H_TH)%gcs     => gcs
    call mult_ini(mult(T_H_TH), Hmatind, myid-1, prim%n_prim, parts)
    mult(T_L_TL)%mult_type = 1
    mult(T_L_TL)%amat    => mat(1:prim%groups_on_node,Trange)
    mult(T_L_TL)%bmat    => mat(1:prim%groups_on_node,Lrange)
    mult(T_L_TL)%cmat    => mat(1:prim%groups_on_node,TLrange)
    mult(T_L_TL)%ahalo   => halo(Trange)
    mult(T_L_TL)%chalo   => halo(TLrange)
    mult(T_L_TL)%ltrans  => ltrans(Trange)
    !mult(T_L_TL)%bindex => Lmatind
    mult(T_L_TL)%parts   => parts
    mult(T_L_TL)%prim    => prim
    mult(T_L_TL)%gcs     => gcs
    call mult_ini(mult(T_L_TL), Lmatind, myid-1, prim%n_prim, parts)
    mult(TL_T_L)%mult_type = 2
    mult(TL_T_L)%amat    => mat(1:prim%groups_on_node,Lrange)
    mult(TL_T_L)%bmat    => mat(1:prim%groups_on_node,TTrrange)
    mult(TL_T_L)%cmat    => mat(1:prim%groups_on_node,TLrange)
    mult(TL_T_L)%ahalo   => halo(Lrange)
    mult(TL_T_L)%chalo   => halo(TLrange)
    mult(TL_T_L)%ltrans  => ltrans(Lrange)
    !mult(TL_T_L)%bindex => TTrmatind
    mult(TL_T_L)%parts   => parts
    mult(TL_T_L)%prim    => prim
    mult(TL_T_L)%gcs     => gcs
    call mult_ini(mult(TL_T_L), TTrmatind, myid-1, prim%n_prim, parts)
    mult(S_X_SX)%mult_type = 1
    mult(S_X_SX)%amat    => mat(1:prim%groups_on_node,Srange)
    mult(S_X_SX)%bmat    => mat(1:prim%groups_on_node,Xrange)
    mult(S_X_SX)%cmat    => mat(1:prim%groups_on_node,SXrange)
    mult(S_X_SX)%ahalo   => halo(Srange)
    mult(S_X_SX)%chalo   => halo(SXrange)
    mult(S_X_SX)%ltrans  => ltrans(Srange)
    mult(S_X_SX)%parts   => parts
    mult(S_X_SX)%prim    => prim
    mult(S_X_SX)%gcs     => gcs
    call mult_ini(mult(S_X_SX), Xmatind, myid-1, prim%n_prim, parts)
    ! 24-34 : for multiplications of atomf-based matrices
    if (atomf.ne.sf) then
       mult(aSa_sCaTr_aSs)%mult_type = 1
       mult(aSa_sCaTr_aSs)%amat    => mat(1:prim%groups_on_node,aSa_range)
       mult(aSa_sCaTr_aSs)%bmat    => mat(1:prim%groups_on_node,SFcoeffTr_range)
       mult(aSa_sCaTr_aSs)%cmat    => mat(1:prim%groups_on_node,aSs_range)
       mult(aSa_sCaTr_aSs)%ahalo   => halo(aSa_range)
       mult(aSa_sCaTr_aSs)%chalo   => halo(aSs_range)
       mult(aSa_sCaTr_aSs)%ltrans  => ltrans(aSa_range)
       !mult(aSa_sCaTr_aSs)%bindex => SFcoeffTr_matind
       mult(aSa_sCaTr_aSs)%parts   => parts
       mult(aSa_sCaTr_aSs)%prim    => prim
       mult(aSa_sCaTr_aSs)%gcs     => gcs
       call mult_ini(mult(aSa_sCaTr_aSs), SFcoeffTr_matind, myid-1, prim%n_prim, parts)

       mult(sCa_aSs_sSs)%mult_type = 1
       mult(sCa_aSs_sSs)%amat    => mat(1:prim%groups_on_node,SFcoeff_range)
       mult(sCa_aSs_sSs)%bmat    => mat(1:prim%groups_on_node,aSs_range)
       mult(sCa_aSs_sSs)%cmat    => mat(1:prim%groups_on_node,Srange)
       mult(sCa_aSs_sSs)%ahalo   => halo(SFcoeff_range)
       mult(sCa_aSs_sSs)%chalo   => halo(Srange)
       mult(sCa_aSs_sSs)%ltrans  => ltrans(SFcoeff_range)
       !mult(sCa_aSs_sSs)%bindex => aSs_matind
       mult(sCa_aSs_sSs)%parts   => parts
       mult(sCa_aSs_sSs)%prim    => prim
       mult(sCa_aSs_sSs)%gcs     => gcs
       call mult_ini(mult(sCa_aSs_sSs), aSs_matind, myid-1, prim%n_prim, parts)

       mult(aHa_sCaTr_aHs)%mult_type = 1
       mult(aHa_sCaTr_aHs)%amat    => mat(1:prim%groups_on_node,aHa_range)
       mult(aHa_sCaTr_aHs)%bmat    => mat(1:prim%groups_on_node,SFcoeffTr_range)
       mult(aHa_sCaTr_aHs)%cmat    => mat(1:prim%groups_on_node,aHs_range)
       mult(aHa_sCaTr_aHs)%ahalo   => halo(aHa_range)
       mult(aHa_sCaTr_aHs)%chalo   => halo(aHs_range)
       mult(aHa_sCaTr_aHs)%ltrans  => ltrans(aHa_range)
       !mult(aHa_sCaTr_aHs)%bindex => SFcoeffTr_matind
       mult(aHa_sCaTr_aHs)%parts   => parts
       mult(aHa_sCaTr_aHs)%prim    => prim
       mult(aHa_sCaTr_aHs)%gcs     => gcs
       call mult_ini(mult(aHa_sCaTr_aHs), SFcoeffTr_matind, myid-1, prim%n_prim, parts)

       mult(sCa_aHs_sHs)%mult_type = 1
       mult(sCa_aHs_sHs)%amat    => mat(1:prim%groups_on_node,SFcoeff_range)
       mult(sCa_aHs_sHs)%bmat    => mat(1:prim%groups_on_node,aHs_range)
       mult(sCa_aHs_sHs)%cmat    => mat(1:prim%groups_on_node,Hrange)
       mult(sCa_aHs_sHs)%ahalo   => halo(SFcoeff_range)
       mult(sCa_aHs_sHs)%chalo   => halo(Hrange)
       mult(sCa_aHs_sHs)%ltrans  => ltrans(SFcoeff_range)
       !mult(sCa_aHs_sHs)%bindex => aHs_matind
       mult(sCa_aHs_sHs)%parts   => parts
       mult(sCa_aHs_sHs)%prim    => prim
       mult(sCa_aHs_sHs)%gcs     => gcs
       call mult_ini(mult(sCa_aHs_sHs), aHs_matind, myid-1, prim%n_prim, parts)

       mult(sCaTr_sSs_aSs)%mult_type = 2
       mult(sCaTr_sSs_aSs)%amat    => mat(1:prim%groups_on_node,aSs_range)
       mult(sCaTr_sSs_aSs)%bmat    => mat(1:prim%groups_on_node,STr_range)
       mult(sCaTr_sSs_aSs)%cmat    => mat(1:prim%groups_on_node,SFcoeffTr_range)
       mult(sCaTr_sSs_aSs)%ahalo   => halo(aSs_range)
       mult(sCaTr_sSs_aSs)%chalo   => halo(SFcoeffTr_range)
       mult(sCaTr_sSs_aSs)%ltrans  => ltrans(aSs_range)
       !mult(sCaTr_sSs_aSs)%bindex => STr_matind
       mult(sCaTr_sSs_aSs)%parts   => parts
       mult(sCaTr_sSs_aSs)%prim    => prim
       mult(sCaTr_sSs_aSs)%gcs     => gcs
       call mult_ini(mult(sCaTr_sSs_aSs), STr_matind, myid-1, prim%n_prim, parts)

       mult(aSs_sCa_aSa)%mult_type = 2
       mult(aSs_sCa_aSa)%amat    => mat(1:prim%groups_on_node,aSa_range)
       mult(aSs_sCa_aSa)%bmat    => mat(1:prim%groups_on_node,SFcoeffTr_range)
       mult(aSs_sCa_aSa)%cmat    => mat(1:prim%groups_on_node,aSs_range)
       mult(aSs_sCa_aSa)%ahalo   => halo(aSa_range)
       mult(aSs_sCa_aSa)%chalo   => halo(aSs_range)
       mult(aSs_sCa_aSa)%ltrans  => ltrans(aSa_range)
       !mult(aSs_sCa_aSa)%bindex => SFcoeffTr_matind
       mult(aSs_sCa_aSa)%parts   => parts
       mult(aSs_sCa_aSa)%prim    => prim
       mult(aSs_sCa_aSa)%gcs     => gcs
       call mult_ini(mult(aSs_sCa_aSa), SFcoeffTr_matind, myid-1, prim%n_prim, parts)

       mult(sCaTr_sHs_aHs)%mult_type = 2
       mult(sCaTr_sHs_aHs)%amat    => mat(1:prim%groups_on_node,aHs_range)
       mult(sCaTr_sHs_aHs)%bmat    => mat(1:prim%groups_on_node,HTr_range)
       mult(sCaTr_sHs_aHs)%cmat    => mat(1:prim%groups_on_node,SFcoeffTr_range)
       mult(sCaTr_sHs_aHs)%ahalo   => halo(aHs_range)
       mult(sCaTr_sHs_aHs)%chalo   => halo(SFcoeffTr_range)
       mult(sCaTr_sHs_aHs)%ltrans  => ltrans(aHs_range)
       !mult(sCaTr_sHs_aHs)%bindex => HTr_matind
       mult(sCaTr_sHs_aHs)%parts   => parts
       mult(sCaTr_sHs_aHs)%prim    => prim
       mult(sCaTr_sHs_aHs)%gcs     => gcs
       call mult_ini(mult(sCaTr_sHs_aHs), HTr_matind, myid-1, prim%n_prim, parts)

       mult(aHs_sCa_aHa)%mult_type = 2
       mult(aHs_sCa_aHa)%amat    => mat(1:prim%groups_on_node,aHa_range)
       mult(aHs_sCa_aHa)%bmat    => mat(1:prim%groups_on_node,SFcoeffTr_range)
       mult(aHs_sCa_aHa)%cmat    => mat(1:prim%groups_on_node,aHs_range)
       mult(aHs_sCa_aHa)%ahalo   => halo(aHa_range)
       mult(aHs_sCa_aHa)%chalo   => halo(aHs_range)
       mult(aHs_sCa_aHa)%ltrans  => ltrans(aHa_range)
       !mult(aHs_sCa_aHa)%bindex => SFcoeffTr_matind
       mult(aHs_sCa_aHa)%parts   => parts
       mult(aHs_sCa_aHa)%prim    => prim
       mult(aHs_sCa_aHa)%gcs     => gcs
       call mult_ini(mult(aHs_sCa_aHa), SFcoeffTr_matind, myid-1, prim%n_prim, parts)

       mult(sSs_sSa_sCa)%mult_type = 2
       mult(sSs_sSa_sCa)%amat    => mat(1:prim%groups_on_node,SFcoeff_range)
       mult(sSs_sSa_sCa)%bmat    => mat(1:prim%groups_on_node,aSs_range)
       mult(sSs_sSa_sCa)%cmat    => mat(1:prim%groups_on_node,Srange)
       mult(sSs_sSa_sCa)%ahalo   => halo(SFcoeff_range)
       mult(sSs_sSa_sCa)%chalo   => halo(Srange)
       mult(sSs_sSa_sCa)%ltrans  => ltrans(SFcoeff_range)
       !mult(sSs_sSa_sCa)%bindex => aSs_matind
       mult(sSs_sSa_sCa)%parts   => parts
       mult(sSs_sSa_sCa)%prim    => prim
       mult(sSs_sSa_sCa)%gcs     => gcs
       call mult_ini(mult(sSs_sSa_sCa), aSs_matind, myid-1, prim%n_prim, parts)

       mult(sHs_sHa_sCa)%mult_type = 2
       mult(sHs_sHa_sCa)%amat    => mat(1:prim%groups_on_node,SFcoeff_range)
       mult(sHs_sHa_sCa)%bmat    => mat(1:prim%groups_on_node,aHs_range)
       mult(sHs_sHa_sCa)%cmat    => mat(1:prim%groups_on_node,Hrange)
       mult(sHs_sHa_sCa)%ahalo   => halo(SFcoeff_range)
       mult(sHs_sHa_sCa)%chalo   => halo(Hrange)
       mult(sHs_sHa_sCa)%ltrans  => ltrans(SFcoeff_range)
       !mult(sHs_sHa_sCa)%bindex => aHs_matind
       mult(sHs_sHa_sCa)%parts   => parts
       mult(sHs_sHa_sCa)%prim    => prim
       mult(sHs_sHa_sCa)%gcs     => gcs
       call mult_ini(mult(sHs_sHa_sCa), aHs_matind, myid-1, prim%n_prim, parts)

       if (flag_LFD) then
          mult(aLa_aSa_aLSa)%mult_type = 1
          mult(aLa_aSa_aLSa)%amat    => mat(1:prim%groups_on_node,LD_range)
          mult(aLa_aSa_aLSa)%bmat    => mat(1:prim%groups_on_node,aSa_range)
          mult(aLa_aSa_aLSa)%cmat    => mat(1:prim%groups_on_node,LD_range)
          mult(aLa_aSa_aLSa)%ahalo   => halo(LD_range)
          mult(aLa_aSa_aLSa)%chalo   => halo(LD_range)
          mult(aLa_aSa_aLSa)%ltrans  => ltrans(LD_range)
          !mult(aLa_aSa_aLSa)%bindex => aSa_matind
          mult(aLa_aSa_aLSa)%parts   => parts
          mult(aLa_aSa_aLSa)%prim    => prim
          mult(aLa_aSa_aLSa)%gcs     => gcs
          call mult_ini(mult(aLa_aSa_aLSa), aSa_matind, myid-1, prim%n_prim, parts)

          mult(aLa_aHa_aLHa)%mult_type = 1
          mult(aLa_aHa_aLHa)%amat    => mat(1:prim%groups_on_node,LD_range)
          mult(aLa_aHa_aLHa)%bmat    => mat(1:prim%groups_on_node,aHa_range)
          mult(aLa_aHa_aLHa)%cmat    => mat(1:prim%groups_on_node,LD_range)
          mult(aLa_aHa_aLHa)%ahalo   => halo(LD_range)
          mult(aLa_aHa_aLHa)%chalo   => halo(LD_range)
          mult(aLa_aHa_aLHa)%ltrans  => ltrans(LD_range)
          !mult(aLa_aHa_aLHa)%bindex => aHa_matind
          mult(aLa_aHa_aLHa)%parts   => parts
          mult(aLa_aHa_aLHa)%prim    => prim
          mult(aLa_aHa_aLHa)%gcs     => gcs
          call mult_ini(mult(aLa_aHa_aLHa), aHa_matind, myid-1, prim%n_prim, parts)
       endif
    endif ! atomf
    if(flag_neutral_atom_projector) then
       ra = rcut(aNArange)
       rc = rcut(aHa_range)
       !if (rc >= ra) then
          mult(aNA_NAa_aHa)%mult_type = 1
          mult(aNA_NAa_aHa)%amat    => mat(1:prim%groups_on_node,aNArange)
          mult(aNA_NAa_aHa)%bmat    => mat(1:prim%groups_on_node,NAarange)
          mult(aNA_NAa_aHa)%cmat    => mat(1:prim%groups_on_node,aHa_range)
          mult(aNA_NAa_aHa)%ahalo   => halo(aNArange)
          mult(aNA_NAa_aHa)%chalo   => halo(aHa_range)
          mult(aNA_NAa_aHa)%ltrans  => ltrans(aNArange)
          !mult(aNA_NAa_aHa)%bindex => NAamatind
          mult(aNA_NAa_aHa)%parts   => parts
          mult(aNA_NAa_aHa)%prim    => prim
          mult(aNA_NAa_aHa)%gcs     => gcs
          call mult_ini(mult(aNA_NAa_aHa), NAamatind, myid-1, prim%n_prim, parts)
       !else
       !   mult(aNA_NAa_aHa)%mult_type = 2
       !   mult(aNA_NAa_aHa)%amat    => mat(1:prim%groups_on_node,aHa_range)
       !   mult(aNA_NAa_aHa)%bmat    => mat(1:prim%groups_on_node,aNArange)
       !   mult(aNA_NAa_aHa)%cmat    => mat(1:prim%groups_on_node,aNArange)
       !   mult(aNA_NAa_aHa)%ahalo   => halo(aHa_range)
       !   mult(aNA_NAa_aHa)%chalo   => halo(aNArange)
       !   mult(aNA_NAa_aHa)%ltrans  => ltrans(aHa_range)
       !   !mult(aNA_NAa_aHa)%bindex => aNAmatind
       !   mult(aNA_NAa_aHa)%parts   => parts
       !   mult(aNA_NAa_aHa)%prim    => prim
       !   mult(aNA_NAa_aHa)%gcs     => gcs
       !   call mult_ini(mult(aNA_NAa_aHa), aNAmatind, myid-1, prim%n_prim, parts)
       !end if
       mult(aHa_aNA_aNA)%mult_type = 2
       mult(aHa_aNA_aNA)%amat    => mat(1:prim%groups_on_node,aNArange)
       mult(aHa_aNA_aNA)%bmat    => mat(1:prim%groups_on_node,NAarange)
       mult(aHa_aNA_aNA)%cmat    => mat(1:prim%groups_on_node,aHa_range)
       mult(aHa_aNA_aNA)%ahalo   => halo(aNArange)
       mult(aHa_aNA_aNA)%chalo   => halo(aHa_range)
       mult(aHa_aNA_aNA)%ltrans  => ltrans(aNArange)
       !mult(aHa_aNA_aNA)%bindex => NAamatind
       mult(aHa_aNA_aNA)%parts   => parts
       mult(aHa_aNA_aNA)%prim    => prim
       mult(aHa_aNA_aNA)%gcs     => gcs
       call mult_ini(mult(aHa_aNA_aNA), NAamatind, myid-1, prim%n_prim, parts)
    end if
    !end if
    ! call stop_timer(tmr_std_matrices)

    return
  end subroutine immi
  !!***

  !!****f* mult_module/fmmi *
  !!
  !!  NAME
  !!   fmmi - finish matrix multiplication indexing
  !!  USAGE
  !!   fmmi
  !!  PURPOSE
  !!   Deallocates all memory associated with matrix multiplication
  !!   indexing
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   03/07/2001
  !!  MODIFICATION HISTORY
  !!   10/08/2001 dave
  !!    Added nullify calls - probably a little fussy, but can't harm
  !!   14:57, 16/11/2004 dave
  !!    Removed TS_TS_T and S_TS_T - not needed
  !!   2008/05/22 ast
  !!    Added timers
  !!   2016/08/09 21:30 nakata
  !!    Added lines for atomf-based matrices and LFD matrices
  !!   2017/03/08 15:00 nakata
  !!    Removed PAOP_PS_H, dSrange, dHrange and PAOPrange which are no longer used
  !!   2017/12/05 10:33 dave with TM and NW (Mizuho)
  !!    Deallocating NA projector matrices
  !!   2018/07/11 12:12 dave
  !!    Adding dissociate_matrices and deallocation of mat for empty bundle fix
  !!  SOURCE
  !!
  subroutine fmmi(prim)

    use basic_types
    use global_module, only: iprint_mat, flag_LFD
    use matrix_module, only: deallocate_comms_data
    use matrix_data
    use maxima_module, only: maxnabaprocs
    use pseudopotential_common, ONLY: flag_neutral_atom_projector

    implicit none

    ! Passed variables
    type(primary_set) :: prim

    ! Local variables
    integer :: i,j

!    call start_timer(tmr_std_matrices)
   ! Multiplications
    nullify(mult(aHa_AP_AP)%amat,mult(aHa_AP_AP)%bmat,mult(aHa_AP_AP)%cmat,       &
    mult(aHa_AP_AP)%ahalo,mult(aHa_AP_AP)%chalo,mult(aHa_AP_AP)%ltrans,           &
    mult(aHa_AP_AP)%bindex,mult(aHa_AP_AP)%parts,mult(aHa_AP_AP)%prim,            &
    mult(aHa_AP_AP)%gcs)
    nullify(mult(L_S_LS)%amat,mult(L_S_LS)%bmat,mult(L_S_LS)%cmat,          &
    mult(L_S_LS)%ahalo,mult(L_S_LS)%chalo,mult(L_S_LS)%ltrans,              &
    mult(L_S_LS)%bindex,mult(L_S_LS)%parts,mult(L_S_LS)%prim,               &
    mult(L_S_LS)%gcs)
    nullify(mult(L_H_LH)%amat,mult(L_H_LH)%bmat,mult(L_H_LH)%cmat,          &
    mult(L_H_LH)%ahalo,mult(L_H_LH)%chalo,mult(L_H_LH)%ltrans,              &
    mult(L_H_LH)%bindex,mult(L_H_LH)%parts,mult(L_H_LH)%prim,               &
    mult(L_H_LH)%gcs)
    nullify(mult(T_S_TS)%amat,mult(T_S_TS)%bmat,mult(T_S_TS)%cmat,          &
    mult(T_S_TS)%ahalo,mult(T_S_TS)%chalo,mult(T_S_TS)%ltrans,              &
    mult(T_S_TS)%bindex,mult(T_S_TS)%parts,mult(T_S_TS)%prim,               &
    mult(T_S_TS)%gcs)
    !nullify(mult(S_TS_S)%amat,mult(S_TS_S)%bmat,mult(S_TS_S)%cmat,         &
    !mult(S_TS_S)%ahalo,mult(S_TS_S)%chalo,mult(S_TS_S)%ltrans,             &
    !mult(S_TS_S)%bindex,mult(S_TS_S)%parts,mult(S_TS_S)%prim,              &
    !mult(S_TS_S)%gcs)
    !nullify(mult(S_LS_L)%amat,mult(S_LS_L)%bmat,mult(S_LS_L)%cmat,         &
    !mult(S_LS_L)%ahalo,mult(S_LS_L)%chalo,mult(S_LS_L)%ltrans,             &
    !mult(S_LS_L)%bindex,mult(S_LS_L)%parts,mult(S_LS_L)%prim,              &
    !mult(S_LS_L)%gcs)
    nullify(mult(SL_S_SLS)%amat,mult(SL_S_SLS)%bmat,mult(SL_S_SLS)%cmat,    &
    mult(SL_S_SLS)%ahalo,mult(SL_S_SLS)%chalo,mult(SL_S_SLS)%ltrans,        &
    mult(SL_S_SLS)%bindex,mult(SL_S_SLS)%parts,mult(SL_S_SLS)%prim,         &
    mult(SL_S_SLS)%gcs)
    nullify(mult(LH_L_SLS)%amat,mult(LH_L_SLS)%bmat,mult(LH_L_SLS)%cmat,    &
    mult(LH_L_SLS)%ahalo,mult(LH_L_SLS)%chalo,mult(LH_L_SLS)%ltrans,        &
    mult(LH_L_SLS)%bindex,mult(LH_L_SLS)%parts,mult(LH_L_SLS)%prim,         &
    mult(LH_L_SLS)%gcs)
    nullify(mult(LS_L_LSL)%amat,mult(LS_L_LSL)%bmat,mult(LS_L_LSL)%cmat,    &
    mult(LS_L_LSL)%ahalo,mult(LS_L_LSL)%chalo,mult(LS_L_LSL)%ltrans,        &
    mult(LS_L_LSL)%bindex,mult(LS_L_LSL)%parts,mult(LS_L_LSL)%prim,         &
    mult(LS_L_LSL)%gcs)
    nullify(mult(HL_S_LSL)%amat,mult(HL_S_LSL)%bmat,mult(HL_S_LSL)%cmat,    &
    mult(HL_S_LSL)%ahalo,mult(HL_S_LSL)%chalo,mult(HL_S_LSL)%ltrans,        &
    mult(HL_S_LSL)%bindex,mult(HL_S_LSL)%parts,mult(HL_S_LSL)%prim,         &
    mult(HL_S_LSL)%gcs)
    nullify(mult(HLS_LS_L)%amat,mult(HLS_LS_L)%bmat,mult(HLS_LS_L)%cmat,    &
    mult(HLS_LS_L)%ahalo,mult(HLS_LS_L)%chalo,mult(HLS_LS_L)%ltrans,        &
    mult(HLS_LS_L)%bindex,mult(HLS_LS_L)%parts,mult(HLS_LS_L)%prim,         &
    mult(HLS_LS_L)%gcs)
    !nullify(mult(LSL_SL_L)%amat,mult(LSL_SL_L)%bmat,mult(LSL_SL_L)%cmat,   &
    !mult(LSL_SL_L)%ahalo,mult(LSL_SL_L)%chalo,mult(LSL_SL_L)%ltrans,       &
    !mult(LSL_SL_L)%bindex,mult(LSL_SL_L)%parts,mult(LSL_SL_L)%prim,        &
    !mult(LSL_SL_L)%gcs)
    nullify(mult(SLS_LS_L)%amat,mult(SLS_LS_L)%bmat,mult(SLS_LS_L)%cmat,    &
    mult(SLS_LS_L)%ahalo,mult(SLS_LS_L)%chalo,mult(SLS_LS_L)%ltrans,        &
    mult(SLS_LS_L)%bindex,mult(SLS_LS_L)%parts,mult(SLS_LS_L)%prim,         &
    mult(SLS_LS_L)%gcs)
    nullify(mult(LHL_SL_S)%amat,mult(LHL_SL_S)%bmat,mult(LHL_SL_S)%cmat,    &
    mult(LHL_SL_S)%ahalo,mult(LHL_SL_S)%chalo,mult(LHL_SL_S)%ltrans,        &
    mult(LHL_SL_S)%bindex,mult(LHL_SL_S)%parts,mult(LHL_SL_S)%prim,         &
    mult(LHL_SL_S)%gcs)
    nullify(mult(LSL_SL_H)%amat,mult(LSL_SL_H)%bmat,mult(LSL_SL_H)%cmat,    &
    mult(LSL_SL_H)%ahalo,mult(LSL_SL_H)%chalo,mult(LSL_SL_H)%ltrans,        &
    mult(LSL_SL_H)%bindex,mult(LSL_SL_H)%parts,mult(LSL_SL_H)%prim,         &
    mult(LSL_SL_H)%gcs)
    nullify(mult(AP_PA_aHa)%amat,mult(AP_PA_aHa)%bmat,mult(AP_PA_aHa)%cmat,       &
    mult(AP_PA_aHa)%ahalo,mult(AP_PA_aHa)%chalo,mult(AP_PA_aHa)%ltrans,           &
    mult(AP_PA_aHa)%bindex,mult(AP_PA_aHa)%parts,mult(AP_PA_aHa)%prim,            &
    mult(AP_PA_aHa)%gcs)
    nullify(mult(TS_T_T)%amat,mult(TS_T_T)%bmat,mult(TS_T_T)%cmat,          &
    mult(TS_T_T)%ahalo,mult(TS_T_T)%chalo,mult(TS_T_T)%ltrans,              &
    mult(TS_T_T)%bindex,mult(TS_T_T)%parts,mult(TS_T_T)%prim,               &
    mult(TS_T_T)%gcs)
    nullify(mult(TH_T_L)%amat,mult(TH_T_L)%bmat,mult(TH_T_L)%cmat,          &
    mult(TH_T_L)%ahalo,mult(TH_T_L)%chalo,mult(TH_T_L)%ltrans,              &
    mult(TH_T_L)%bindex,mult(TH_T_L)%parts,mult(TH_T_L)%prim,               &
    mult(TH_T_L)%gcs)
    nullify(mult(T_H_TH)%amat,mult(T_H_TH)%bmat,mult(T_H_TH)%cmat,          &
    mult(T_H_TH)%ahalo,mult(T_H_TH)%chalo,mult(T_H_TH)%ltrans,              &
    mult(T_H_TH)%bindex,mult(T_H_TH)%parts,mult(T_H_TH)%prim,               &
    mult(T_H_TH)%gcs)
    nullify(mult(T_L_TL)%amat,mult(T_L_TL)%bmat,mult(T_L_TL)%cmat,          &
    mult(T_L_TL)%ahalo,mult(T_L_TL)%chalo,mult(T_L_TL)%ltrans,              &
    mult(T_L_TL)%bindex,mult(T_L_TL)%parts,mult(T_L_TL)%prim,               &
    mult(T_L_TL)%gcs)
    nullify(mult(TL_T_L)%amat,mult(TL_T_L)%bmat,mult(TL_T_L)%cmat,          &
    mult(TL_T_L)%ahalo,mult(TL_T_L)%chalo,mult(TL_T_L)%ltrans,              &
    mult(TL_T_L)%bindex,mult(TL_T_L)%parts,mult(TL_T_L)%prim,               &
    mult(TL_T_L)%gcs)
    nullify(mult(S_X_SX)%amat,mult(S_X_SX)%bmat,mult(S_X_SX)%cmat,          &
    mult(S_X_SX)%ahalo,mult(S_X_SX)%chalo,mult(S_X_SX)%ltrans,              &
    mult(S_X_SX)%bindex,mult(S_X_SX)%parts,mult(S_X_SX)%prim,               &
    mult(S_X_SX)%gcs)
    if (atomf.ne.sf) then
       nullify(mult(aSa_sCaTr_aSs)%amat,mult(aSa_sCaTr_aSs)%bmat,mult(aSa_sCaTr_aSs)%cmat, &
       mult(aSa_sCaTr_aSs)%ahalo,mult(aSa_sCaTr_aSs)%chalo,mult(aSa_sCaTr_aSs)%ltrans,     &
       mult(aSa_sCaTr_aSs)%bindex,mult(aSa_sCaTr_aSs)%parts,mult(aSa_sCaTr_aSs)%prim,      &
       mult(aSa_sCaTr_aSs)%gcs)
       nullify(mult(sCa_aSs_sSs)%amat,mult(sCa_aSs_sSs)%bmat,mult(sCa_aSs_sSs)%cmat, &
       mult(sCa_aSs_sSs)%ahalo,mult(sCa_aSs_sSs)%chalo,mult(sCa_aSs_sSs)%ltrans,     &
       mult(sCa_aSs_sSs)%bindex,mult(sCa_aSs_sSs)%parts,mult(sCa_aSs_sSs)%prim,      &
       mult(sCa_aSs_sSs)%gcs)
       nullify(mult(aHa_sCaTr_aHs)%amat,mult(aHa_sCaTr_aHs)%bmat,mult(aHa_sCaTr_aHs)%cmat, &
       mult(aHa_sCaTr_aHs)%ahalo,mult(aHa_sCaTr_aHs)%chalo,mult(aHa_sCaTr_aHs)%ltrans,     &
       mult(aHa_sCaTr_aHs)%bindex,mult(aHa_sCaTr_aHs)%parts,mult(aHa_sCaTr_aHs)%prim,      &
       mult(aHa_sCaTr_aHs)%gcs)
       nullify(mult(sCa_aHs_sHs)%amat,mult(sCa_aHs_sHs)%bmat,mult(sCa_aHs_sHs)%cmat, &
       mult(sCa_aHs_sHs)%ahalo,mult(sCa_aHs_sHs)%chalo,mult(sCa_aHs_sHs)%ltrans,     &
       mult(sCa_aHs_sHs)%bindex,mult(sCa_aHs_sHs)%parts,mult(sCa_aHs_sHs)%prim,      &
       mult(sCa_aHs_sHs)%gcs)
       nullify(mult(sCaTr_sSs_aSs)%amat,mult(sCaTr_sSs_aSs)%bmat,mult(sCaTr_sSs_aSs)%cmat, &
       mult(sCaTr_sSs_aSs)%ahalo,mult(sCaTr_sSs_aSs)%chalo,mult(sCaTr_sSs_aSs)%ltrans,     &
       mult(sCaTr_sSs_aSs)%bindex,mult(sCaTr_sSs_aSs)%parts,mult(sCaTr_sSs_aSs)%prim,      &
       mult(sCaTr_sSs_aSs)%gcs)
       nullify(mult(aSs_sCa_aSa)%amat,mult(aSs_sCa_aSa)%bmat,mult(aSs_sCa_aSa)%cmat, &
       mult(aSs_sCa_aSa)%ahalo,mult(aSs_sCa_aSa)%chalo,mult(aSs_sCa_aSa)%ltrans,     &
       mult(aSs_sCa_aSa)%bindex,mult(aSs_sCa_aSa)%parts,mult(aSs_sCa_aSa)%prim,      &
       mult(aSs_sCa_aSa)%gcs)
       nullify(mult(sCaTr_sHs_aHs)%amat,mult(sCaTr_sHs_aHs)%bmat,mult(sCaTr_sHs_aHs)%cmat, &
       mult(sCaTr_sHs_aHs)%ahalo,mult(sCaTr_sHs_aHs)%chalo,mult(sCaTr_sHs_aHs)%ltrans,     &
       mult(sCaTr_sHs_aHs)%bindex,mult(sCaTr_sHs_aHs)%parts,mult(sCaTr_sHs_aHs)%prim,      &
       mult(sCaTr_sHs_aHs)%gcs)
       nullify(mult(aHs_sCa_aHa)%amat,mult(aHs_sCa_aHa)%bmat,mult(aHs_sCa_aHa)%cmat, &
       mult(aHs_sCa_aHa)%ahalo,mult(aHs_sCa_aHa)%chalo,mult(aHs_sCa_aHa)%ltrans,     &
       mult(aHs_sCa_aHa)%bindex,mult(aHs_sCa_aHa)%parts,mult(aHs_sCa_aHa)%prim,      &
       mult(aHs_sCa_aHa)%gcs)
       nullify(mult(sSs_sSa_sCa)%amat,mult(sSs_sSa_sCa)%bmat,mult(sSs_sSa_sCa)%cmat, &
       mult(sSs_sSa_sCa)%ahalo,mult(sSs_sSa_sCa)%chalo,mult(sSs_sSa_sCa)%ltrans,     &
       mult(sSs_sSa_sCa)%bindex,mult(sSs_sSa_sCa)%parts,mult(sSs_sSa_sCa)%prim,      &
       mult(sSs_sSa_sCa)%gcs)
       nullify(mult(sHs_sHa_sCa)%amat,mult(sHs_sHa_sCa)%bmat,mult(sHs_sHa_sCa)%cmat, &
       mult(sHs_sHa_sCa)%ahalo,mult(sHs_sHa_sCa)%chalo,mult(sHs_sHa_sCa)%ltrans,     &
       mult(sHs_sHa_sCa)%bindex,mult(sHs_sHa_sCa)%parts,mult(sHs_sHa_sCa)%prim,      &
       mult(sHs_sHa_sCa)%gcs)
       if (flag_LFD) then
          nullify(mult(aLa_aSa_aLSa)%amat,mult(aLa_aSa_aLSa)%bmat,mult(aLa_aSa_aLSa)%cmat, &
          mult(aLa_aSa_aLSa)%ahalo,mult(aLa_aSa_aLSa)%chalo,mult(aLa_aSa_aLSa)%ltrans,     &
          mult(aLa_aSa_aLSa)%bindex,mult(aLa_aSa_aLSa)%parts,mult(aLa_aSa_aLSa)%prim,      &
          mult(aLa_aSa_aLSa)%gcs)
          nullify(mult(aLa_aHa_aLHa)%amat,mult(aLa_aHa_aLHa)%bmat,mult(aLa_aHa_aLHa)%cmat, &
          mult(aLa_aHa_aLHa)%ahalo,mult(aLa_aHa_aLHa)%chalo,mult(aLa_aHa_aLHa)%ltrans,     &
          mult(aLa_aHa_aLHa)%bindex,mult(aLa_aHa_aLHa)%parts,mult(aLa_aHa_aLHa)%prim,      &
          mult(aLa_aHa_aLHa)%gcs)
       endif
    endif ! atomf
    if( flag_neutral_atom_projector ) then
       nullify(mult(aNA_NAa_aHa)%amat,mult(aNA_NAa_aHa)%bmat,mult(aNA_NAa_aHa)%cmat,       &
            mult(aNA_NAa_aHa)%ahalo,mult(aNA_NAa_aHa)%chalo,mult(aNA_NAa_aHa)%ltrans,           &
            mult(aNA_NAa_aHa)%bindex,mult(aNA_NAa_aHa)%parts,mult(aNA_NAa_aHa)%prim,            &
            mult(aNA_NAa_aHa)%gcs)
       nullify(mult(aHa_aNA_aNA)%amat,mult(aHa_aNA_aNA)%bmat,mult(aHa_aNA_aNA)%cmat,       &
            mult(aHa_aNA_aNA)%ahalo,mult(aHa_aNA_aNA)%chalo,mult(aHa_aNA_aNA)%ltrans,           &
            mult(aHa_aNA_aNA)%bindex,mult(aHa_aNA_aNA)%parts,mult(aHa_aNA_aNA)%prim,            &
            mult(aHa_aNA_aNA)%gcs)
    end if
    ! Comms
    call deallocate_comms_data(mult(aHa_AP_AP)%comms)
    call deallocate_comms_data(mult(L_S_LS)%comms)
    call deallocate_comms_data(mult(L_H_LH)%comms)
    call deallocate_comms_data(mult(T_S_TS)%comms)
    !call deallocate_comms_data(mult(S_TS_S)%comms)
    !call deallocate_comms_data(mult(S_LS_L)%comms)
    call deallocate_comms_data(mult(SL_S_SLS)%comms)
    call deallocate_comms_data(mult(LH_L_SLS)%comms)
    call deallocate_comms_data(mult(LS_L_LSL)%comms)
    call deallocate_comms_data(mult(HL_S_LSL)%comms)
    call deallocate_comms_data(mult(HLS_LS_L)%comms)
    !call deallocate_comms_data(mult(LSL_SL_L)%comms)
    call deallocate_comms_data(mult(SLS_LS_L)%comms)
    call deallocate_comms_data(mult(LHL_SL_S)%comms)
    call deallocate_comms_data(mult(LSL_SL_H)%comms)
    call deallocate_comms_data(mult(AP_PA_aHa)%comms)
    call deallocate_comms_data(mult(TS_T_T)%comms)
    call deallocate_comms_data(mult(TH_T_L)%comms)
    call deallocate_comms_data(mult(T_H_TH)%comms)
    call deallocate_comms_data(mult(T_L_TL)%comms)
    call deallocate_comms_data(mult(TL_T_L)%comms)
    if (atomf.ne.sf) then
       call deallocate_comms_data(mult(aSa_sCaTr_aSs)%comms)
       call deallocate_comms_data(mult(sCa_aSs_sSs)%comms)
       call deallocate_comms_data(mult(aHa_sCaTr_aHs)%comms)
       call deallocate_comms_data(mult(sCa_aHs_sHs)%comms)
       call deallocate_comms_data(mult(sCaTr_sSs_aSs)%comms)
       call deallocate_comms_data(mult(aSs_sCa_aSa)%comms)
       call deallocate_comms_data(mult(sCaTr_sHs_aHs)%comms)
       call deallocate_comms_data(mult(aHs_sCa_aHa)%comms)
       call deallocate_comms_data(mult(sSs_sSa_sCa)%comms)
       call deallocate_comms_data(mult(sHs_sHa_sCa)%comms)
       if (flag_LFD) then
          call deallocate_comms_data(mult(aLa_aSa_aLSa)%comms)
          call deallocate_comms_data(mult(aLa_aHa_aLHa)%comms)
       endif
    endif
    if( flag_neutral_atom_projector ) then
       call deallocate_comms_data(mult(aNA_NAa_aHa)%comms)
       call deallocate_comms_data(mult(aHa_aNA_aNA)%comms)
    end if
    do i=1,mx_trans
       do j = 1,gtrans(i)%n_rem_node!maxnabaprocs+1
          !if(associated(pairs(j,i)%submat)) deallocate(pairs(j,i)%submat)
          deallocate(pairs(j,i)%submat)
       end do
    end do
    deallocate(pairs)
    deallocate(Spairind, Lpairind, APpairind, LSpairind, LHpairind, &
               LSLpairind, Tpairind)
    if (atomf.ne.sf) deallocate(aSs_pairind, aHs_pairind, SFcoeff_pairind)
    if(flag_neutral_atom_projector) then
       deallocate(aNApairind)
       if(atomf.ne.sf) deallocate(aNAapairind)
    end if
    call dissociate_matrices
    ! Matrices
    call end_ops(prim,Srange,Smatind,S_trans)
    call end_ops(prim,Lrange,Lmatind,L_trans)
    call end_ops(prim,LTrrange,LTrmatind)
    call end_ops(prim,Hrange,Hmatind)
    call end_ops(prim,APrange, APmatind,AP_trans)
    call end_ops(prim,PArange,PAmatind)
    call end_ops(prim,LSrange, LSmatind,LS_trans)
    call end_ops(prim,SLrange,SLmatind)
    call end_ops(prim,LHrange, LHmatind,LH_trans)
    call end_ops(prim,HLrange,HLmatind)
    call end_ops(prim,LSLrange,LSLmatind,LSL_trans)
    call end_ops(prim,SLSrange,SLSmatind)
    call end_ops(prim,Trange,Tmatind,T_trans)
    call end_ops(prim,TTrrange,TTrmatind)
    call end_ops(prim,TSrange,TSmatind)
    call end_ops(prim,THrange,THmatind)
    call end_ops(prim,TLrange,TLmatind)
    if (atomf.ne.sf) then
       call end_ops(prim,aSa_range,aSa_matind)
       call end_ops(prim,aHa_range,aHa_matind)
       call end_ops(prim,STr_range,STr_matind)
       call end_ops(prim,HTr_range,HTr_matind)
       call end_ops(prim,aSs_range,aSs_matind,aSs_trans)
       call end_ops(prim,aHs_range,aHs_matind,aHs_trans)
       call end_ops(prim,sSa_range,sSa_matind)
       call end_ops(prim,sHa_range,sHa_matind)
       call end_ops(prim,SFcoeff_range,SFcoeff_matind,SFcoeff_trans)
       call end_ops(prim,SFcoeffTr_range,SFcoeffTr_matind)
       if (flag_LFD) call end_ops(prim,LD_range,LD_matind)
    endif
    if( flag_neutral_atom_projector ) then
       call end_ops(prim,aNArange, aNAmatind,aNA_trans)
       call end_ops(prim,NAarange,NAamatind)
    end if
!    call stop_timer(tmr_std_matrices)
    deallocate(mat)
    return
  end subroutine fmmi
  !!***

  !!****f* mult_module/LNV_matrix_multiply *
  !!
  !!  NAME
  !!   main_matrix_multiply
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Performs most matrix multiplications
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   CMG/EHH/IJB/D.R.Bowler
  !!  CREATION DATE
  !!   Unknown
  !!  MODIFICATION HISTORY
  !!   08/06/2001 dave
  !!    Added ROBODoc header and changed to use GenComms
  !!   30/05/2002 dave
  !!    Made all long-range (LH or more) matrices allocatable
  !!   10:43, 06/03/2003 drb
  !!    Removed alias of prim to bundle (ifc problem)
  !!   16:28, 10/05/2005 dave
  !!    Added timing calls and a LOCAL iprint_mult to activate them: we
  !!    can make this global later
  !!   2006/10/19 16:04 dave
  !!    Changed iprint_mult to iprint_mat (global variable)
  !!   Thursday, 2011/07/28  L.Tong
  !!    Adding spin dependent calculations, adding an optional parameter
  !!    "spin" to control which spin channel is to be calculated
  !!    If "spin" parameter is present then the subroutine assumes spin
  !!    polarised calculation and only calculates the relevant
  !!    quantities for the particular spin channel indicated by "spin".
  !!    If "spin" is not present then the subroutine assumes spin
  !!    non-polarised calculation and calculates the relevant quantities
  !!    (with exception to M12) for the TOTAL electron density, which
  !!    has factor of two.
  !!   Wednesday, 2011/08/03 L.Tong
  !!    Fixed the bug caused by the conditional if (present (spin)
  !!    .and. spin == 2), on some compilers which does not short circuit
  !!    and evaluates both sides of .and. and when the subroutine is
  !!    called without the optional spin parameter, then spin == 2 gives
  !!    a segmentation fault. Fixed this by introducing logical variable
  !!    do_spin_down, we work out this variable correctly first, and
  !!    then apply conditional on it in place of the previous present
  !!    (spin) .and. spin == 2 conditional.
  !!   2012/03/04 L.Tong
  !!   - Completely changed implementation of spin. Now spin is regarded
  !!     as an array of dimension nspin. This subroutine the computes
  !!     the necessary matrices for all spin components.
  !!   - mat_M12, mat_M3, mat_M4 and mat_phi are now OPTIONAL variables.
  !!   - electrons and energy are arrays of dimension nspin.
  !!   - mat_M3, mat_M4 and mat_phi are all calculated as matrices
  !!     for single spin channel, even for spin non-polarised case.
  !!   - electrons(nspin) and energy(nspin) now also contain the
  !!     associated number of electrons and band energies in each spin
  !!     channels. And for spin non-polarised case, this means
  !!     electrons(1) and energy(1) contain half of the total number.
  !!   - Added optional variable spin, which if present makes the
  !!     subroutine to only calculate the matrices for a single spin
  !!     channel
  !!   2015/06/08 lat
  !!    - Added experimental backtrace
  !!  SOURCE
  !!
  subroutine LNV_matrix_multiply(electrons, energy, doK, doM1, doM2, &
                                 doM3, doM4, dophi, doE, mat_M12,    &
                                 mat_M3, mat_M4, mat_phi, spin, level)

    use numbers
    use datatypes
    use global_module, only: iprint_mat, IPRINT_TIME_THRES3, nspin, &
                             flag_SpinDependentSF
    use matrix_data,   only: LHrange, SLSrange, LSLrange, Lrange,   &
                             Srange, Hrange
    use GenComms,      only: gsum, inode, ionode, my_barrier, gmax, &
                             mtime
    use timer_module

    implicit none

    ! Passed variables
    real(double), dimension(:) :: electrons, energy
    integer,      dimension(:), optional :: mat_M12, mat_M3, mat_M4, mat_phi
    integer,                    optional :: spin
    integer,                    optional :: level
    integer :: myid
    logical :: doK, doM1, doM2, doM3, doM4, dophi, doE
    
    ! Local variables
    real(double)              :: electrons_1, electrons_2, t0, t1
    integer, dimension(nspin) :: matLH, matHL, matLHL, matHLS, matSLH, &
                                 matA, matSLS, matLSL, matLHLSL,       &
                                 matLSLHL, matB, matC, matSLSLS,       &
                                 matLSLSL
    type(cq_timer)            :: tmr_l_tmp1
    type(cq_timer)            :: backtrace_timer
    integer                   :: backtrace_level
    integer                   :: ss, ss_start, ss_end, ss_SF

    ! This routine is basically calls to matrix operations, so these
    ! are timed within the routines themselves

!****lat<$
    if (       present(level) ) backtrace_level = level+1
    if ( .not. present(level) ) backtrace_level = -100
    call start_backtrace(t=backtrace_timer,who='LNV_matrix_multiply',&
         where=area,level=backtrace_level)
!****lat>$

    call start_timer(tmr_l_tmp1, WITH_LEVEL)

    if (present(spin)) then
       ss_start = spin
       ss_end = spin
    else
       ss_start = 1
       ss_end = nspin
    end if

    if ((doM1 .or. doM2) .and. (.not. present(mat_M12))) &
         call cq_abort("LNVmm needs proper matrix for M12")
    if (doM3 .and. (.not. present(mat_M3))) &
         call cq_abort("LNVmm needs proper matrix for M3")
    if (doM4 .and. (.not. present(mat_M4))) &
         call cq_abort("LNVmm needs proper matrix for M4")
    if (dophi .and. (.not. present(mat_phi))) &
         call cq_abort("LNVmm needs proper matrix for phi")
    myid = inode - 1
    if (iprint_mat > 4) t0 = mtime()

    do ss = ss_start, ss_end
       ss_SF = 1
       if (flag_SpinDependentSF) ss_SF = ss
 
       if (iprint_mat > 4 .and. nspin == 2) then
          if (inode == ionode) write (io_lun, '(1x,"For spin = ",i1," :")') ss
       end if
       call matrix_product(matL(ss), matS(ss_SF), matLS(ss), mult(L_S_LS))
       if (iprint_mat > 4) then
          t1 = mtime()
          if (inode == ionode) write (io_lun, *) 'LS time: ', t1 - t0
          t0 = t1
       end if
       call matrix_transpose(matLS(ss), matSL(ss))

       if (iprint_mat > 4) then
          t1 = mtime()
          if (inode == ionode) write (io_lun, *) 'LS trans time: ', t1 - t0
          t0 = t1
       end if
       if (doM1 .OR. doM2 .OR. doM3) then
          matLH(ss) = allocate_temp_matrix(LHrange, LH_trans)
          matHL(ss) = allocate_temp_matrix(LHrange, LH_trans)
          if (iprint_mat > 4) then
             t1 = mtime()
             if (inode == ionode) write (io_lun, *) 'LH alloc time: ', t1 - t0
             t0 = t1
          end if
          call matrix_product(matL(ss), matH(ss), matLH(ss), mult(L_H_LH))
          if (iprint_mat > 4) then
             t1 = mtime()
             if (inode == ionode) write (io_lun, *) 'LH time: ', t1 - t0
             t0 = t1
          end if
          call matrix_transpose(matLH(ss), matHL(ss))
          if (iprint_mat > 4) then
             t1 = mtime()
             if (inode == ionode) write (io_lun, *) 'LH trans time: ', t1 - t0
             t0 = t1
          end if
       end if ! doM1 .OR. doM2 .OR. doM3
       if (doM1 .OR. doM2) then ! This is type 1 or 2 depending on H
          matLHL(ss) = allocate_temp_matrix(SLSrange, 0)
          call matrix_product(matLH(ss), matL(ss), matLHL(ss), &
               mult(LH_L_SLS))
       endif
       if (doM2) then
          matLHLSL(ss) = allocate_temp_matrix(Srange, S_trans)
          matLSLHL(ss) = allocate_temp_matrix(Srange, S_trans)
          call matrix_product(matLHL(ss), matLS(ss), &
               matLHLSL(ss), mult(LHL_SL_S))
          call matrix_transpose(matLHLSL(ss), matLSLHL(ss))
       endif
       if (doM3) then
          if (iprint_mat > 4) t0 = mtime()
          matHLS(ss) = allocate_temp_matrix(LSLrange, LSL_trans)
          matSLH(ss) = allocate_temp_matrix(LSLrange, LSL_trans)
          matA(ss) = allocate_temp_matrix(LSLrange, LSL_trans)
          matB(ss) = allocate_temp_matrix(Lrange, L_trans)
          matC(ss) = allocate_temp_matrix(Lrange, L_trans)
          if (iprint_mat > 4) then
             t1 = mtime()
             if (inode == ionode) write (io_lun, *) 'Alloc time: ', t1 - t0
             t0 = t1
          end if
          call matrix_product(matHL(ss), matS(ss_SF), matHLS(ss), mult(HL_S_LSL))
          if (iprint_mat > 4) then
             t1 = mtime()
             if (inode == ionode) write (io_lun, *) 'HLS time: ', t1 - t0
             t0 = t1
          end if
          call matrix_product(matHL(ss), matS(ss_SF), matHLS(ss), mult(HL_S_LSL))
          if (iprint_mat > 4) then
             t1 = mtime()
             if (inode == ionode) write (io_lun, *) 'HLS time: ', t1 - t0
             t0 = t1
          end if
          call matrix_transpose(matHLS(ss), matSLH(ss))
          if (iprint_mat > 4) then
             t1 = mtime()
             if (inode == ionode) write (io_lun, *) 'HLS trans time: ', t1 - t0
             t0 = t1
          end if
          call matrix_sum(zero, matA(ss), two, matHLS(ss))
          call matrix_sum(one, matA(ss), one, matSLH(ss))
          if (iprint_mat > 4) then
             t1 = mtime()
             if (inode == ionode) write (io_lun, *) 'A time: ', t1 - t0
             t0 = t1
          end if
          call matrix_product(matA(ss), matSL(ss), matB(ss), &
               mult(HLS_LS_L))
          if (iprint_mat > 4) then
             t1 = mtime()
             if (inode == ionode) write (io_lun, *) 'B time: ', t1 - t0
             t0 = t1
          end if
          call matrix_transpose(matB(ss), matC(ss))
          if (iprint_mat > 4) then
             t1 = mtime()
             if (inode == ionode) write (io_lun, *) 'B trans time: ', t1 - t0
             t0 = t1
          end if
       endif
       if (dophi) then
          if (iprint_mat > 4) t0 = mtime()
          matSLS(ss) = allocate_temp_matrix(SLSrange, 0)
          matSLSLS(ss) = allocate_temp_matrix(Lrange, 0)
          if (iprint_mat > 4) then
             t1 = mtime()
             if (inode == ionode) write (io_lun, *) 'SLS alloc time: ', t1 - t0
             t0 = t1
          end if
          call matrix_product(matSL(ss), matS(ss_SF), matSLS(ss), &
                              mult(SL_S_SLS))
          if (iprint_mat > 4) then
             t1 = mtime()
             if (inode == ionode) write (io_lun, *) 'SLS time: ', t1 - t0
             t0 = t1
          end if
          call matrix_product(matSLS(ss), matSL(ss), matSLSLS(ss), &
                              mult(SLS_LS_L))
          if (iprint_mat > 4) then
             t1 = mtime()
             if (inode == ionode) write (io_lun, *) 'SLSLS time: ', t1 - t0
             t0 = t1
          end if
       end if
       if (doK .OR. doE .OR. doM4) then
          if (iprint_mat > 4) t0 = mtime()
          matLSL(ss) = allocate_temp_matrix(LSLrange, 0)
          matLSLSL(ss) = allocate_temp_matrix(Hrange, 0)
          if (iprint_mat > 4) then
             t1 = mtime()
             if (inode == ionode) write (io_lun, *) 'LSL alloc time: ', t1 - t0
             t0 = t1
          end if
          call matrix_product(matLS(ss), matL(ss), matLSL(ss), &
               mult(LS_L_LSL))
          if (iprint_mat > 4) then
             t1 = mtime()
             if (inode == ionode) write (io_lun, *) 'LSL time: ', t1 - t0
             t0 = t1
          end if
          call matrix_product(matLSL(ss), matLS(ss), &
               matLSLSL(ss), mult(LSL_SL_H))
          if (iprint_mat > 4) then
             t1 = mtime()
             if (inode == ionode) write (io_lun, *) 'LSLSL time: ', t1 - t0
             t0 = t1
          end if
       end if
       if (doK .OR. doE) then
          ! work out matK or matK_dn
          call matrix_sum(zero, matK(ss), three, matLSL(ss))
          call matrix_sum(one, matK(ss), -two, matLSLSL(ss))
          if (iprint_mat > 4) then
             t1 = mtime()
             if (inode == ionode) write (io_lun, *) 'K time: ', t1 - t0
             t0 = t1
          end if
       end if
       if (doM4) then
          ! if spin parameter is present, means we are doing spin
          ! polarised calculations, and we do not have the factor
          ! of two used in spin non-polarised case.
          call matrix_sum(zero, mat_M4(ss), 12.0_double, matLSL(ss))
          call matrix_sum(one, mat_M4(ss), -12.0_double, matLSLSL(ss))
       end if
       if (doE) then
          if (iprint_mat > 4) t0 = mtime()
          energy(ss) = matrix_product_trace(matH(ss), matK(ss))
          if (iprint_mat > 4) then
             t1 = mtime()
             if (inode == ionode) write (io_lun, *) 'energy dot time: ', t1 - t0
             t0 = t1
          end if
       end if
       if (doM1 .AND. doM2) then
          call matrix_sum(zero, mat_M12(ss), three, matLHL(ss))
          call matrix_sum(one,  mat_M12(ss), -two,  matLHLSL(ss))
          call matrix_sum(one,  mat_M12(ss), -two,  matLSLHL(ss))
       endif
       if (doM3) then
          if (iprint_mat > 4) t0 = mtime()
          call matrix_sum(zero, mat_M3(ss), three, matHLS(ss))
          if (iprint_mat > 4) then
             t1 = mtime()
             if (inode == ionode) write (io_lun, *) 'M3 time: ', t1 - t0
             t0 = t1
          end if
          call matrix_sum(one, mat_M3(ss), three, matSLH(ss))
          if (iprint_mat > 4) then
             t1 = mtime()
             if (inode == ionode) write (io_lun, *) 'M3 add time: ', t1 - t0
             t0 = t1
          end if
          call matrix_sum(one, mat_M3(ss), -one, matB(ss))
          if (iprint_mat > 4) then
             t1 = mtime()
             if (inode == ionode) write (io_lun, *) 'M3 axpy B time: ', t1 - t0
             t0 = t1
          end if
          call matrix_sum(one, mat_M3(ss), -one, matC(ss))
          if (iprint_mat > 4) then
             t1 = mtime()
             if (inode == ionode) write (io_lun, *) 'M3 axpy C time: ', t1 - t0
             t0 = t1
          end if
       end if
       if (dophi) then
          if (iprint_mat > 4) t0 = mtime()
          call matrix_sum(zero, mat_phi(ss), one, matSLS(ss))
          if (iprint_mat > 4) then
             t1 = mtime()
             if (inode == ionode) write (io_lun, *) 'SLS prune time: ', t1 - t0
             t0 = t1
          end if
          electrons_1 = matrix_product_trace(matL(ss), mat_phi(ss))
          electrons_2 = matrix_product_trace(matL(ss), matSLSLS(ss))
          electrons(ss) = three * electrons_1 - two * electrons_2
          if (iprint_mat > 4) then
             t1 = mtime()
             if (inode == ionode) &
                  write (io_lun, *) 'electrons dot time: ', t1 - t0
             t0 = t1
          end if
          call matrix_sum(six, mat_phi(ss), -six, matSLSLS(ss))
          if (iprint_mat > 4) then
             t1 = mtime()
             if (inode == ionode) write (io_lun, *) 'phi time: ', t1 - t0
             t0 = t1
          end if
       end if ! dophi
       
       ! free temporary matrices
       call my_barrier
       if (iprint_mat > 4) t0 = mtime()
       if (doK .or. doE .or. doM4) then
          call free_temp_matrix(matLSLSL(ss))
          call free_temp_matrix(matLSL(ss))
       end if
       if (dophi) then
          call free_temp_matrix(matSLSLS(ss))
          call free_temp_matrix(matSLS(ss))
       end if
       if (doM3) then
          call free_temp_matrix(matC(ss))
          call free_temp_matrix(matB(ss))
          call free_temp_matrix(matA(ss))
          call free_temp_matrix(matSLH(ss))
          call free_temp_matrix(matHLS(ss))
       end if
       if (doM2) then
          call free_temp_matrix(matLSLHL(ss))
          call free_temp_matrix(matLHLSL(ss))
       end if
       if (doM1 .or. doM2) then
          call free_temp_matrix(matLHL(ss))
       end if
       if (doM1 .or. doM2 .or. doM3) then
          call free_temp_matrix(matHL(ss))
          call free_temp_matrix(matLH(ss))
       end if
       if (iprint_mat > 4) then
          t1 = mtime()
          if (inode == ionode) write (io_lun, *) 'dealloc time: ', t1 - t0
       end if
           
    end do ! ss (spin)
    
    call stop_print_timer(tmr_l_tmp1, "LNV_matrix_multiply", &
                          IPRINT_TIME_THRES3)
    
!****lat<$
    call stop_backtrace(t=backtrace_timer,who='LNV_matrix_multiply')
!****lat>$

    return
  end subroutine LNV_matrix_multiply
  !!***


  !!****f* mult_module/symmetrise_L *
  !!
  !!  NAME
  !!   symmetrise_L
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Ensures that L is symmetric
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   CMG/D.R.Bowler
  !!  CREATION DATE
  !! ?
  !!  MODIFICATION HISTORY
  !!   08/06/2001 dave
  !!    Added ROBODoc header, used BLAS routines
  !!   Tuesday, 2011/08/02 L.Tong
  !!    Added support for spin
  !!   2012/03/04 L.Tong
  !!    Changed spin implementation. Now the optional parameter spin
  !!    controls if individual spin components are to be symmetrised,
  !!    if spin is not present, then all spin channels are
  !!    symmetrised.
  !!  SOURCE
  !!
  subroutine symmetrise_L(spin)

    use datatypes
    use numbers
    use matrix_data,   only: Lrange
    use global_module, only: nspin

    implicit none

    ! Passed variables
    ! spin = 1 for up, spin = 2 for dn, not used = no spin or all spin
    ! dependents on whether spin polarisation is on
    integer, optional :: spin

    ! Local variables
    integer :: length, mattmp, ss

    mattmp = allocate_temp_matrix(Lrange, L_trans)
    if (present(spin)) then
       if (spin == 2) then
          call matrix_transpose(matL(2), mattmp)
          call matrix_sum(half, matL(2), half, mattmp)
       else if (spin == 1) then
          call matrix_transpose(matL(1), mattmp)
          call matrix_sum(half, matL(1), half, mattmp)
       else
          call cq_abort('symmetrise_L: spin is undefined ', spin)
       end if
    else
       do ss = 1, nspin
          call matrix_transpose(matL(ss), mattmp)
          call matrix_sum(half, matL(ss), half, mattmp)
       end do
    end if
    call free_temp_matrix(mattmp)
   
    return
  end subroutine symmetrise_L
  !!***

  !!****f* mult_module/symmetrise_mat *
  !!
  !!  NAME
  !!   symmetrise_mat
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Symmetrises the matrix designated in dummy arguments
  !!  INPUTS
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   Michiaki Arita
  !!  CREATION DATE
  !!   2013/08/23
  !!  MODIFICATION HISTORY
  !!   2019/11/13 tsuyoshi
  !!     integer matA -> integer matA(n_matrix) for spin polarisation
  !!  SOURCE
  !!
  subroutine symmetrise_mat(range,n_matrix,trans,matA)

    ! Module usage
    use datatypes
    use numbers

    implicit none

    ! passed variables
    !integer :: range,trans,matA
    integer :: range,trans,n_matrix
    integer :: matA(n_matrix)

    ! local variables
    integer :: mattmp
    integer :: isize

    mattmp = allocate_temp_matrix(range,trans)
    do isize = 1, n_matrix
     call matrix_transpose(matA(isize),mattmp)
     call matrix_sum(half,matA(isize),half,mattmp)
    enddo
    call free_temp_matrix(mattmp)

    return
  end subroutine symmetrise_mat
  !!***

  !!****f* mult_module/end_ops *
  !!
  !!  NAME
  !!   end_ops
  !!  USAGE
  !!   end_ops(matrix range, global transpose)
  !!  PURPOSE
  !!   Deallocates memory associated with matrices and their indexing
  !!  INPUTS
  !!   integer :: range - matrix to deallocate
  !!   integer, OPTIONAL :: gtran - global transpose to deallocate
  !!  USES
  !!
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   04/07/2001
  !!  MODIFICATION HISTORY
  !!
  !!  SOURCE
  !!
  subroutine end_ops(prim, range, matind, gtran)

    use basic_types
    use matrix_data, only: mat, halo

    implicit none

    ! Passed variables
    type(primary_set) :: prim
    integer           :: range
    integer, pointer, dimension(:) :: matind
    integer, OPTIONAL :: gtran
    ! Local variables
    integer :: i

    ! Global transpose
    if (present(gtran)) call deallocate_trans_rem(gtrans(gtran))
    ! Matrices
    call deallocate_trans(ltrans(range))
    call deallocate_halo(halo(range))
    call deallocate_matrix(mat(:,range),prim%groups_on_node)
    deallocate(matind)
    do i = 1, prim%groups_on_node
       nullify(mat(i,range)%n_nab, mat(i,range)%i_seq,    &
               mat(i,range)%i_acc, mat(i,range)%i_nd_acc, &
               mat(i,range)%npxyz, mat(i,range)%ndimj)
    end do

    return
  end subroutine end_ops
  !!***


  !*** New routines start here ***!
  
  !!****f* mult_module/associate_matrices *
  !!
  !!  NAME
  !!   associate_matrices
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Associates the matrix pointers in private type with actual
  !!   matrices in matrix_data
  !!
  !!   I realise that this is ugly; it does, however, mean that we can
  !!   test the new matrix routines without too many changes
  !!  INPUTS
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   2006/01/25 17:04 dave
  !!  MODIFICATION HISTORY
  !!   2008/05/22 ast
  !!    Added timers
  !!   2011/03/28 L.Tong
  !!    Added assignment indices for matrices: matH_dn, matL_dn,
  !!    matK_dn, matdH_dn, matSL_dn, matL_dnS and matphi_dn
  !!   2011/11/24 L.Tong
  !!    Added assignment indices for matrices:
  !!    matU_dn and matUT_dn
  !!   2012/03/04 L.Tong
  !!    Changed spin implementation, now mat(1) corresponds to spin up
  !!    and mat(2) corresponds to spin down. The global variable nspin
  !!    = 1 | 2 controls if the calculation will be spin polarised or
  !!    not
  !!   2016/08/09 21:30 nakata
  !!    Added lines for atomf-based matrices and LFD matrices
  !!   2017/03/08 15:00 nakata
  !!    Removed dSrange, dHrange, matdS and matdH which are no longer used
  !!   2017/12/05 10:34 dave
  !!    Associating NA-projector matrices
  !!   2018/11/13 17:30 nakata
  !!    Changed matS, matT, matTtran, matKE, matNL and matNA to be spin_SF dependent
  !!  SOURCE
  !!
  subroutine associate_matrices

    use numbers,       only: zero
    use matrix_data,   only: Srange, Hrange, Lrange, Trange, PArange, &
                             LSrange, APrange, SLrange, TTrrange,     &
                             Xrange, SXrange, aSa_range, aHa_range,   &
                             SFcoeff_range, SFcoeffTr_range,          &
                             mat, aNArange, NAarange
    use global_module, only: area_matrices, iprint_mat, nspin, nspin_SF
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use GenComms,      only: inode, ionode
!%%!    use matrix_data, only: data_S, data_H, data_L, data_T, data_Ttran, &
!%%!         data_LS, data_SL, data_SC, data_CS, data_K, data_M4, data_U, data_phi, &
!%%!         data_M12, data_UT, data_KE, data_NL, Srange, Hrange, Lrange, Trange, PArange,&
!%%!         LSrange, APrange, SLrange, dSrange, dHrange, data_dH, data_dS, mat,TTrrange
    use pseudopotential_common, ONLY: flag_neutral_atom_projector

    implicit none

    integer :: stat, i, spin, spin_SF
    logical :: allocated_tags = .false.

    ! Allocate spin dependent matrices
    if (.not. allocated_tags) then
       allocate(matH(nspin), matL(nspin),   matLS(nspin),  matSL(nspin), &
                matK(nspin), matphi(nspin), matM12(nspin), matM4(nspin), &
                matU(nspin), matUT(nspin),  matX(nspin),   matSX(nspin), STAT=stat)
       if (stat /= 0) &
            call cq_abort('associate_matrices: failed to allocate spin &
                           &depdendent matrix tags', nspin, stat)
       ! For atomf-based calculations
       allocate(matHatomf(nspin), matKatomf(nspin), matXatomf(nspin), STAT=stat)
       if (stat /= 0) &
            call cq_abort('associate_matrices: failed to allocate spin &
                           &depdendent ATOMF matrix tags', nspin, stat)
       ! For spin_SF dependent matrices
       if (atomf.ne.sf) then
          allocate(matSFcoeff(nspin_SF),  matSFcoeff_tran(nspin_SF), &
                   matdSFcoeff(nspin_SF), matdSFcoeff_e(nspin_SF), STAT=stat)
          if (stat /= 0) &
               call cq_abort('associate_matrices: failed to allocate spin_SF &
                              &depdendent SFcoefficient matrix tags', nspin_SF, stat)
       endif
       allocate(matS(nspin_SF), matT(nspin_SF), matTtran(nspin_SF), &
                matKE(nspin_SF), matNL(nspin_SF), matNA(nspin_SF), STAT=stat) 
       if (stat /= 0) &
            call cq_abort('associate_matrices: failed to allocate spin_SF &
                           &depdendent matrix tags', nspin_SF, stat)
       allocated_tags = .true.
    end if

    ! Assign indices for the matrices that we need
    matS(1)     = 1
    matH(1)     = 2
    matL(1)     = 3
    matT(1)     = 4
    matTtran(1) = 5
    matLS(1)    = 6
    matSL(1)    = 7
    matAP       = 8
    matPA       = 9
    matK(1)     = 10
    matM4(1)    = 11
    matU(1)     = 12
    matphi(1)   = 13
    matM12(1)   = 14
    matUT(1)    = 15
    matKE(1)    = 16
    matNL(1)    = 17
    matX(1)     = 18
    matSX(1)    = 19
    current_matrix = 19
    if (nspin == 2) then
       matH(2)   = 20
       matL(2)   = 21
       matK(2)   = 22
       matLS(2)  = 23
       matSL(2)  = 24
       matphi(2) = 25
       matU(2)   = 26
       matUT(2)  = 27
       matM12(2) = 28
       matM4(2)  = 29
       matX(2)   = 30
       matSX(2)  = 31
       current_matrix = 31
    endif
    if (nspin_SF == 2) then
       ! when support functions are spin-dependent
       matS(2)     = current_matrix + 1
       matT(2)     = current_matrix + 2
       matTtran(2) = current_matrix + 3
       matKE(2)    = current_matrix + 4
       matNL(2)    = current_matrix + 5
       current_matrix = current_matrix + 5
    endif
    if (atomf.eq.sf) then
    ! when atomic functions are euivalent to SFs (blips and one_to_one PAOs)
       matSatomf    = matS(1)
       matKEatomf   = matKE(1)
       matNLatomf   = matNL(1)
       matHatomf(1) = matH(1)
       matKatomf(1) = matK(1)
       matXatomf(1) = matX(1)
       if (nspin == 2) then
          matHatomf(2) = matH(2)
          matKatomf(2) = matK(2)
          matXatomf(2) = matX(2)
       endif
    else
    ! contracted SFs
       matSatomf          = current_matrix + 1     ! 20 or 32
       matHatomf(1)       = current_matrix + 2     ! 21 or 33
       matKatomf(1)       = current_matrix + 3     ! 22 or 34
       matKEatomf         = current_matrix + 4     ! 23 or 35
       matNLatomf         = current_matrix + 5     ! 24 or 36
       matXatomf(1)       = current_matrix + 6     ! 25 or 37
       matSFcoeff(1)      = current_matrix + 7     ! 26 or 38
       matSFcoeff_tran(1) = current_matrix + 8     ! 27 or 39
       matdSFcoeff(1)     = current_matrix + 9     ! 28 or 40
       matdSFcoeff_e(1)   = current_matrix + 10    ! 29 or 41
       current_matrix     = current_matrix + 10    ! 29 or 41
       if (nspin == 2) then
          matHatomf(2)       = current_matrix + 1  ! 42
          matKatomf(2)       = current_matrix + 2  ! 43
          matXatomf(2)       = current_matrix + 3  ! 44
          current_matrix     = current_matrix + 3  ! 44
          if (nspin_SF == 2) then
             matSFcoeff(2)      = current_matrix + 1  ! 45
             matSFcoeff_tran(2) = current_matrix + 2  ! 46
             matdSFcoeff(2)     = current_matrix + 3  ! 47
             matdSFcoeff_e(2)   = current_matrix + 4  ! 48
             current_matrix     = current_matrix + 4  ! 48
          endif
       endif
    endif ! atomf
    if( flag_neutral_atom_projector ) then
       current_matrix = current_matrix + 1
       mataNA = current_matrix
       current_matrix = current_matrix + 1
       matNAa = current_matrix
       current_matrix = current_matrix + 1
       matNA(1)  = current_matrix
       if (nspin_SF == 2) then
          current_matrix = current_matrix + 1
          matNA(2)  = current_matrix
       endif
       if (atomf.eq.sf) then
          matNAatomf = matNA(1)
       else
          current_matrix = current_matrix + 1
          matNAatomf  = current_matrix
       end if
    end if

!%%!! Now associate pointers with correct arrays
!%%!    mat_p(matS  )%matrix => data_S
!%%!    mat_p(matH  )%matrix => data_H
!%%!    mat_p(matL  )%matrix => data_L
!%%!    mat_p(matT  )%matrix => data_T
!%%!    mat_p(matTtran)%matrix => data_Ttran
!%%!    mat_p(matLS )%matrix => data_LS
!%%!    mat_p(matSL )%matrix => data_SL
!%%!    mat_p(matAP )%matrix => data_AP
!%%!    mat_p(matPA )%matrix => data_PA
!%%!    mat_p(matK  )%matrix => data_K
!%%!    mat_p(matM4 )%matrix => data_M4
!%%!    mat_p(matU  )%matrix => data_U
!%%!    mat_p(matphi)%matrix => data_phi
!%%!    mat_p(matM12)%matrix => data_M12
!%%!    mat_p(matUT )%matrix => data_UT
!%%!    mat_p(matKE )%matrix => data_KE
!%%!    mat_p(matNL )%matrix => data_NL

    ! Set the dimensions of the arrays
    ! Type for f1
    mat_p(matAP   )%sf1_type = atomf
    mat_p(matPA   )%sf1_type = nlpf
    do spin_SF = 1, nspin_SF
       mat_p(matS(spin_SF)    )%sf1_type = sf
       mat_p(matT(spin_SF)    )%sf1_type = sf
       mat_p(matTtran(spin_SF))%sf1_type = sf    
       mat_p(matKE(spin_SF)   )%sf1_type = sf
       mat_p(matNL(spin_SF)   )%sf1_type = sf
    enddo
    do spin = 1, nspin
       mat_p(matL(spin)  )%sf1_type = sf
       mat_p(matK(spin)  )%sf1_type = sf
       mat_p(matH(spin)  )%sf1_type = sf
       mat_p(matphi(spin))%sf1_type = sf
       mat_p(matM12(spin))%sf1_type = sf
       mat_p(matM4(spin) )%sf1_type = sf
       mat_p(matU(spin)  )%sf1_type = atomf
       mat_p(matUT(spin) )%sf1_type = nlpf
       mat_p(matLS(spin) )%sf1_type = sf    
       mat_p(matSL(spin) )%sf1_type = sf
       mat_p(matX(spin)  )%sf1_type = sf
       mat_p(matSX(spin) )%sf1_type = sf
    end do
    if (atomf.ne.sf) then
       mat_p(matSatomf )%sf1_type = atomf
       mat_p(matKEatomf)%sf1_type = atomf
       mat_p(matNLatomf)%sf1_type = atomf
       do spin = 1, nspin
          mat_p(matKatomf(spin))%sf1_type = atomf
          mat_p(matHatomf(spin))%sf1_type = atomf
          mat_p(matXatomf(spin))%sf1_type = atomf
       enddo
       do spin_SF = 1, nspin_SF
          mat_p(matSFcoeff(spin_SF)     )%sf1_type = sf
          mat_p(matSFcoeff_tran(spin_SF))%sf1_type = atomf
          mat_p(matdSFcoeff(spin_SF)    )%sf1_type = sf
          mat_p(matdSFcoeff_e(spin_SF)  )%sf1_type = sf
       enddo
    endif
    if( flag_neutral_atom_projector ) then
       mat_p(mataNA  )%sf1_type = atomf
       mat_p(matNAa  )%sf1_type = napf
       do spin_SF = 1, nspin_SF
          mat_p(matNA(spin_SF))%sf1_type = sf
       enddo
       if (atomf.ne.sf) then
          mat_p(matNAatomf   )%sf1_type = atomf
       end if
    end if
    ! Type for f2
    mat_p(matAP   )%sf2_type = nlpf
    mat_p(matPA   )%sf2_type = atomf
    do spin_SF = 1, nspin_SF
       mat_p(matS(spin_SF)    )%sf2_type = sf
       mat_p(matT(spin_SF)    )%sf2_type = sf
       mat_p(matTtran(spin_SF))%sf2_type = sf
       mat_p(matKE(spin_SF)   )%sf2_type = sf
       mat_p(matNL(spin_SF)   )%sf2_type = sf
    enddo
    do spin = 1, nspin
       mat_p(matL(spin)  )%sf2_type = sf
       mat_p(matK(spin)  )%sf2_type = sf
       mat_p(matH(spin)  )%sf2_type = sf
       mat_p(matphi(spin))%sf2_type = sf
       mat_p(matM12(spin))%sf2_type = sf
       mat_p(matM4(spin) )%sf2_type = sf
       mat_p(matU(spin)  )%sf2_type = nlpf
       mat_p(matUT(spin) )%sf2_type = atomf
       mat_p(matLS(spin) )%sf2_type = sf
       mat_p(matSL(spin) )%sf2_type = sf
       mat_p(matX(spin)  )%sf2_type = sf
       mat_p(matSX(spin) )%sf2_type = sf
    end do
    if (atomf.ne.sf) then
       mat_p(matSatomf )%sf2_type = atomf
       mat_p(matKEatomf)%sf2_type = atomf
       mat_p(matNLatomf)%sf2_type = atomf
       do spin = 1, nspin
          mat_p(matKatomf(spin))%sf2_type = atomf
          mat_p(matHatomf(spin))%sf2_type = atomf
          mat_p(matXatomf(spin))%sf2_type = atomf
       enddo
       do spin_SF = 1, nspin_SF
          mat_p(matSFcoeff(spin_SF)     )%sf2_type = atomf
          mat_p(matSFcoeff_tran(spin_SF))%sf2_type = sf
          mat_p(matdSFcoeff(spin_SF)    )%sf2_type = atomf
          mat_p(matdSFcoeff_e(spin_SF)  )%sf2_type = atomf
       enddo
    endif
    if( flag_neutral_atom_projector ) then
       mat_p(mataNA  )%sf2_type = napf
       mat_p(matNAa  )%sf2_type = atomf
       do spin_SF = 1, nspin_SF
          mat_p(matNA(spin_SF))%sf2_type = atomf
       enddo
       if (atomf.ne.sf) then
          mat_p(matNAatomf   )%sf2_type = atomf
       end if
    end if
    ! Set up index translation for ranges
    matrix_index(matAP   ) = APrange
    matrix_index(matPA   ) = PArange
    do spin_SF = 1, nspin_SF
       matrix_index(matS(spin_SF)    ) = Srange
       matrix_index(matT(spin_SF)    ) = Trange
       matrix_index(matTtran(spin_SF)) = TTrrange
       matrix_index(matKE(spin_SF)   ) = Hrange
       matrix_index(matNL(spin_SF)   ) = Hrange    
    enddo
    do spin = 1, nspin
       matrix_index(matL(spin)  ) = Lrange
       matrix_index(matK(spin)  ) = Hrange
       matrix_index(matH(spin)  ) = Hrange
       matrix_index(matphi(spin)) = Lrange
       matrix_index(matM12(spin)) = Srange
       matrix_index(matM4(spin) ) = Srange
       matrix_index(matU(spin)  ) = APrange
       matrix_index(matUT(spin) ) = PArange
       matrix_index(matLS(spin) ) = LSrange
       matrix_index(matSL(spin) ) = SLrange
       matrix_index(matX(spin)  ) = Hrange  ! Xrange
       matrix_index(matSX(spin) ) = SXrange
    end do
    if (atomf.ne.sf) then
       matrix_index(matSatomf ) = aSa_range
       matrix_index(matKEatomf) = aHa_range
       matrix_index(matNLatomf) = aHa_range
       do spin = 1, nspin
          matrix_index(matKatomf(spin)) = aHa_range
          matrix_index(matHatomf(spin)) = aHa_range
          matrix_index(matXatomf(spin)) = aHa_range   ! Xatomf_range
       enddo
       do spin_SF = 1, nspin_SF
          matrix_index(matSFcoeff(spin_SF)     ) = SFcoeff_range
          matrix_index(matSFcoeff_tran(spin_SF)) = SFcoeffTr_range
          matrix_index(matdSFcoeff(spin_SF)    ) = SFcoeff_range
          matrix_index(matdSFcoeff_e(spin_SF)  ) = SFcoeff_range
       end do
    endif
    if( flag_neutral_atom_projector ) then
       matrix_index(mataNA  ) = aNArange
       matrix_index(matNAa  ) = NAarange
       do spin_SF = 1, nspin_SF
          matrix_index(matNA(spin_SF)) = Hrange    
       enddo
       if (atomf.ne.sf) then
          matrix_index(matNAatomf) = aHa_range
       end if
    end if
    ! Real lengths
    ! Length
!%%!    mat_p(matS  )%length = nsf*nsf*mx_at_prim*mx_snab
!%%!    mat_p(matH  )%length = nsf*nsf*mx_at_prim*mx_hnab
!%%!    mat_p(matL  )%length = nsf*nsf*mx_at_prim*mx_lnab
!%%!    mat_p(matT  )%length = nsf*nsf*mx_at_prim*mx_tnab
!%%!    mat_p(matTtran)%length = nsf*nsf*mx_at_prim*mx_tnab
!%%!    mat_p(matLS )%length = nsf*nsf*mx_at_prim*mx_lsnab
!%%!    mat_p(matSL )%length = nsf*nsf*mx_at_prim*mx_lsnab
!%%!    mat_p(matAP )%length = natomf*nlpf*mx_at_prim*mx_scnab
!%%!    mat_p(matPA )%length = nlpf*natomf*mx_at_prim*mx_scnab
!%%!    mat_p(matK  )%length = nsf*nsf*mx_at_prim*mx_hnab
!%%!    mat_p(matM4 )%length = nsf*nsf*mx_at_prim*mx_snab
!%%!    mat_p(matU  )%length = nsf*nlpf*mx_at_prim*mx_scnab
!%%!    mat_p(matphi)%length = nsf*nsf*mx_at_prim*mx_lnab
!%%!    mat_p(matM12)%length = nsf*nsf*mx_at_prim*mx_snab
!%%!    mat_p(matUT )%length = nlpf*nsf*mx_at_prim*mx_scnab
!%%!    mat_p(matKE )%length = nsf*nsf*mx_at_prim*mx_hnab
!%%!    mat_p(matNL )%length = nsf*nsf*mx_at_prim*mx_hnab

    ! Allocalting the core matrices
    do i = 1, current_matrix
       if (mat_p(i)%length == 0) then
          mat_p(i)%length = mat(1,matrix_index(i))%length
          call start_timer(tmr_std_allocation)
          allocate(mat_p(i)%matrix(mat_p(i)%length),STAT=stat)
          if(stat/=0) call cq_abort("Error allocating matrix ",i,mat_p(i)%length)
          call reg_alloc_mem(area_matrices,mat_p(i)%length,type_dbl)
          call stop_timer(tmr_std_allocation)
          mat_p(i)%matrix = zero
       !else if(mat(1,matrix_index(i))%length>mat_p(i)%length) then
       else !if(mat(1,matrix_index(i))%length/=mat_p(i)%length) then
          call start_timer(tmr_std_allocation)
          deallocate(mat_p(i)%matrix,STAT=stat)
          if(stat/=0) call cq_abort("Error deallocating matrix ",i,mat_p(i)%length)
          nullify(mat_p(i)%matrix)
          call reg_dealloc_mem(area_matrices,mat_p(i)%length,type_dbl)
          call stop_timer(tmr_std_allocation)
          mat_p(i)%length = mat(1,matrix_index(i))%length
          call start_timer(tmr_std_allocation)
          allocate(mat_p(i)%matrix(mat_p(i)%length),STAT=stat)
          if(stat/=0) call cq_abort("Error allocating matrix ",i,mat_p(i)%length)
          call reg_alloc_mem(area_matrices,mat_p(i)%length,type_dbl)
          call stop_timer(tmr_std_allocation)
          mat_p(i)%matrix = zero
       end if
    end do
    if(iprint_mat > 3 .AND. inode == ionode) then
       do stat=1,current_matrix
          write (io_lun,fmt='(8x,"Proc: ",i5," Matrix no., length: ",i24)')&
                stat, mat_p(stat)%length
       end do
    end if
    ! Set up index translation for transposes
    trans_index(matAP   ) = AP_trans
    trans_index(matPA   ) = AP_trans
    do spin_SF = 1, nspin_SF
       trans_index(matS(spin_SF)    ) = S_trans
       trans_index(matT(spin_SF)    ) = T_trans
       trans_index(matTtran(spin_SF)) = T_trans
       trans_index(matKE(spin_SF)   ) = 0
       trans_index(matNL(spin_SF)   ) = 0
    enddo
    do spin = 1, nspin
       trans_index(matL(spin)  ) = L_trans
       trans_index(matK(spin)  ) = 0
       trans_index(matH(spin)  ) = 0
       trans_index(matphi(spin)) = L_trans
       trans_index(matM12(spin)) = S_trans
       trans_index(matM4(spin) ) = S_trans
       trans_index(matU(spin)  ) = AP_trans
       trans_index(matUT(spin) ) = AP_trans
       trans_index(matLS(spin) ) = LS_trans
       trans_index(matSL(spin) ) = LS_trans
       trans_index(matX(spin)  ) = 0 ! S_trans
       trans_index(matSX(spin) ) = 0
    end do
    if (atomf.ne.sf) then
       trans_index(matSatomf ) = 0
       trans_index(matKEatomf) = 0
       trans_index(matNLatomf) = 0
       do spin = 1, nspin
          trans_index(matKatomf(spin)) = 0
          trans_index(matHatomf(spin)) = 0 ! Hatomf_trans
          trans_index(matXatomf(spin)) = 0 ! Satomf_trans
       enddo
       do spin_SF = 1, nspin_SF
          trans_index(matSFcoeff(spin_SF)     ) = SFcoeff_trans 
          trans_index(matSFcoeff_tran(spin_SF)) = SFcoeff_trans 
          trans_index(matdSFcoeff(spin_SF)    ) = SFcoeff_trans 
          trans_index(matdSFcoeff_e(spin_SF)  ) = SFcoeff_trans 
       end do
    endif
    if( flag_neutral_atom_projector ) then
       trans_index(mataNA  ) = aNA_trans
       trans_index(matNAa  ) = aNA_trans
       do spin_SF = 1, nspin_SF
          trans_index(matNA(spin_SF))  = S_trans
       enddo
       if (atomf.ne.sf) then
          trans_index(matNAatomf) = aNAa_trans
       end if
    end if

    return
  end subroutine associate_matrices
  !!***


  !!****f* mult_module/dissociate_matrices
  !! PURPOSE
  !!   Deallocate memories
  !! INPUTS
  !! OUTPUT
  !! RETURN VALUE
  !! AUTHOR
  !!   drb
  !! CREATION DATE
  !!
  !! MODIFICATION HISTORY
  !!   2011/07/19 L.Tong
  !!     Added lines for spin polarisation
  !!   2012/03/28 L.Tong
  !!   - Changed spin implementation
  !!   2016/08/09 21:30 nakata
  !!    Added lines for atomf-based matrices and LFD matrices
  !!   2017/03/08 15:00 nakata
  !!    Removed matdS and matdH which are no longer used
  !!   2017/12/05 10:41 dave
  !!    Changed to use loop (as associate_matrices does for allocation)
  !! SOURCE
  !!
  subroutine dissociate_matrices

    use global_module, only: area_matrices
    use memory_module, only: reg_dealloc_mem, type_dbl

    implicit none

    integer :: stat, i

    call start_timer(tmr_std_allocation)
    do i=current_matrix,1,-1
       deallocate(mat_p(i)%matrix,STAT=stat) 
       if (stat /= 0) call cq_abort("Error deallocating matrix ",i,mat_p(i)%length)
       nullify(mat_p(i)%matrix)
       call reg_dealloc_mem(area_matrices,mat_p(i)%length,type_dbl)
       mat_p(i)%length = 0
    end do
    call stop_timer(tmr_std_allocation)

    return
  end subroutine dissociate_matrices
  !!*****


  !!****f* mult_module/mult_wrap *
  !!
  !!  NAME
  !!   mult_wrap
  !!  USAGE
  !!
  !!  PURPOSE
  !!   sbrt mult_wrap: Acts as a wrapper for matrix multiply to reorder the
  !!   arguments - so we can make the same call for strong and weak
  !!   extension or reduction - either way (provided B or B^T is passed
  !!   correctly) the result of the operation is C=A.B.  N.B. trans depends
  !!   on mult type - if 1, Atrans, if 2, Ctrans
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   26/04/00
  !!  MODIFICATION HISTORY
  !!   08/06/2001 dave
  !!    Added ROBODoc header
  !!   2008/05/22 ast
  !!    Added timers
  !!  SOURCE
  !!
  subroutine matrix_product(A, B, C, mult,debug)

    ! Module usage
    use matrix_module,   only: matrix_mult
    use multiply_module, only: mat_mult
    use GenComms,        only: myid

    implicit none

    ! Passed variables
    integer :: A, B, C
    integer, OPTIONAL :: debug
    type(matrix_mult) :: mult

!    call start_timer(tmr_std_matrices)
    ! Type one does 3=1.2 - result to 3
    ! Type two does 1=3.2^T
    !write(io_lun,*) 'Length of A: ',size(mat_p(A)%matrix),mat_p(A)%length,A
    !write(io_lun,*) 'Length of B: ',size(mat_p(B)%matrix),mat_p(B)%length,B
    !write(io_lun,*) 'Length of C: ',size(mat_p(C)%matrix),mat_p(C)%length,C
    !write(io_lun,*) 'Matrices: ',A,B,C,matrix_index(A),matrix_index(B),matrix_index(C)
    if (present(debug)) then
       !write(io_lun,*) 'Debug set in matrix_product'
       if(mult%mult_type==1) then
          call mat_mult(myid, mat_p(A)%matrix, mat_p(A)%length, &
                        mat_p(B)%matrix, mat_p(B)%length, &
                        mat_p(C)%matrix, mat_p(C)%length, &
                        mult, debug)
       else
          call mat_mult(myid, mat_p(C)%matrix, mat_p(C)%length, &
                        mat_p(B)%matrix, mat_p(B)%length, &
                        mat_p(A)%matrix, mat_p(A)%length, &
                        mult, debug)
       endif
     else
        if (mult%mult_type == 1) then
          call mat_mult(myid, mat_p(A)%matrix, mat_p(A)%length, &
                        mat_p(B)%matrix, mat_p(B)%length, &
                        mat_p(C)%matrix, mat_p(C)%length, &
                        mult)
       else
          call mat_mult(myid, mat_p(C)%matrix, mat_p(C)%length, &
                        mat_p(B)%matrix, mat_p(B)%length, &
                        mat_p(A)%matrix, mat_p(A)%length, &
                        mult)
       end if
    end if
    ! call stop_timer(tmr_std_matrices)

    return
  end subroutine matrix_product
  !!***

  subroutine matrix_transpose(A, AT)

    use mpi
    use GenComms,        only: myid, cq_abort
    use matrix_data,     only: halo
    use global_module,   only: sf
    use maxima_module,   only: maxatomsproc, maxnabaprocs
    use multiply_module, only: loc_trans, mat_tran
    use comms_module,    only: send_trans_data

    implicit none

    ! Passed variables
    integer :: A, AT
    ! Local variables
    integer :: mat_temp
    integer :: isend, ii, ierr=0, stat
    integer, allocatable, dimension(:)  :: nreqs
    integer, dimension(MPI_STATUS_SIZE) :: mpi_stat

!    call start_timer(tmr_std_matrices)
    ! Allocate number of requests
    allocate(nreqs(maxnabaprocs+1), STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error allocating number of maximum requests &
                        &in matrix_transpose: ", maxnabaprocs, stat)
    ! Allocate temporary matrix to receive local transpose
    if (mat_p(A)%sf2_type /= sf .or. mat_p(A)%sf1_type /= sf) then
       mat_temp = allocate_temp_matrix(matrix_index(A), 0, &
                                       mat_p(A)%sf2_type,  &
                                       mat_p(A)%sf1_type)
    else
       mat_temp = allocate_temp_matrix(matrix_index(A),0)
    end if
    ! Do local transpose for A
    !write(io_lun,*) myid,' Doing loc_trans'
    call loc_trans(ltrans(matrix_index(A)), halo(matrix_index(A)), &
                   mat_p(A)%matrix, mat_p(A)%length,               &
                   mat_p(mat_temp)%matrix, mat_p(mat_temp)%length, 0)
    ! Send out the data
    !write(io_lun,*) myid,' Doing send_trans'
    call send_trans_data(maxatomsproc, gtrans(trans_index(A)),         &
                         mat_p(mat_temp)%matrix, myid+1, isend, nreqs, &
                         pairs(:,trans_index(A)), maxnabaprocs+1)
    ! Receive it and built AT
    !write(io_lun,*) myid,' Doing mat_tran'
    call mat_tran(maxatomsproc, myid+1, gtrans(trans_index(A)), &
                  mat_p(AT)%matrix, mat_p(mat_temp)%matrix,     &
                  pairs(:,trans_index(A)), maxnabaprocs+1)

    !MPI_Wait  T. Miyazaki 2003/11/18
    do ii= 1, isend
       call MPI_Wait(nreqs(ii), mpi_stat, ierr)
       if (ierr /= 0) &
            call cq_abort("Error in mat_trans: waiting for send to finish", ii)
    end do
    call free_temp_matrix(mat_temp)
    ! Free the array of requests
    deallocate(nreqs, STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error deallocating requests in matrix_transpose: ", stat)
    ! call stop_timer(tmr_std_matrices)

    return    
  end subroutine matrix_transpose


  !! This is equivalent to axpy: A = alpha*A + beta*B
  subroutine matrix_sum(alpha, A, beta, B)

    use datatypes
    use multiply_module, only: matrix_add
    use matrix_data,     only: mat
    use primary_module,  only: bundle
    use GenComms,        only: myid
    use GenBlas,         only: axpy, scal

    implicit none

    ! Passed variables
    integer      :: A, B
    real(double) :: alpha, beta

!    call start_timer(tmr_std_matrices)
    !write(io_lun,*) 'Summing matrices, ranges: ',A,B,matrix_index(A),matrix_index(B)
    if (matrix_index(A) == matrix_index(B)) then ! The matrices are the same range
       !write(io_lun,*) myid,' Scaling ',size(mat_p(A)%matrix),mat_p(A)%length
       call scal(mat_p(A)%length, alpha, mat_p(A)%matrix, 1)
       !mat_p(A)%matrix = alpha*mat_p(A)%matrix
       !write(io_lun,*) myid,' Adding ',size(mat_p(A)%matrix),size(mat_p(B)%matrix),mat_p(A)%length,mat_p(B)%length
       call axpy(mat_p(A)%length, beta, mat_p(B)%matrix, 1, &
                 mat_p(A)%matrix, 1)
    else
       if (matrix_index(A) < 1 .or. matrix_index(A) > mx_matrices) &
            call cq_abort("Matrix error in matrix_sum: ", &
                          matrix_index(A), A)
       if (matrix_index(B) < 1 .or. matrix_index(B) > mx_matrices) &
            call cq_abort("Matrix error in matrix_sum: ", &
                          matrix_index(B), B)
       call scal(mat_p(A)%length, alpha, mat_p(A)%matrix, 1)
       call matrix_add(alpha, mat_p(A)%matrix, mat_p(A)%length,       &
                       mat(:,matrix_index(A)), beta, mat_p(B)%matrix, &
                       mat_p(B)%length, mat(:,matrix_index(B)),       &
                       bundle, myid)
    end if
!    call stop_timer(tmr_std_matrices)
    return
  end subroutine matrix_sum


  subroutine matrix_scale(alpha, A)

    use datatypes
    use GenBlas, only: scal

    implicit none

    integer :: A
    real(double) :: alpha

    ! call start_timer(tmr_std_matrices)
    ! We could do this with a call to scal (the BLAS routine)
    if (A < 1 .or. A > current_matrix) &
         call cq_abort("Error in matrix_scale: ", A, mx_matrices)
    !write(io_lun,*) 'Scaling: ',A,alpha,size(mat_p(A)%matrix)
    call scal(mat_p(A)%length, alpha, mat_p(A)%matrix, 1)
    !mat_p(A)%matrix = alpha*mat_p(A)%matrix
    ! call stop_timer(tmr_std_matrices)

    return
  end subroutine matrix_scale


  !!****f* mult_module/Trace *
  !!
  !!  NAME
  !!   Trace
  !!  USAGE
  !!   Trace(matrix)
  !!  PURPOSE
  !!   Calculates the trace of the specified matrix
  !!
  !!   Note: in writing this, we have assumed that we'll only ever
  !!   take the trace of a matrix which is <sf|A|sf> (i.e. that it has
  !!   support functions on-site).  This is explicit in the use of
  !!   nsf_species for the inner loop.  It would be perfectly possible
  !!   to change this with a switch based on the type of matrix if
  !!   that became necessary
  !!  INPUTS
  !!   integer :: A (matrix index)
  !!  USES
  !!
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!
  !!  MODIFICATION HISTORY
  !!   22/06/2001 dave
  !!    Changed MPI_allreduce to gsum
  !!   2008/05/22 ast
  !!    Added timers
  !!   2024/05/27 nakata
  !!    bugfix: corrected sf to sf1 for loc to check overflow
  !!  SOURCE
  !!
  real(double) function matrix_trace(A)

    use datatypes
    use numbers
    use matrix_data,    only: mat
    use primary_module, only: bundle
    use GenComms,       only: gsum

    implicit none

    ! Passed variables
    integer :: A

    ! Local variables
    integer :: np, i, j, iprim, sf1, loc, Ah

!    call start_timer(tmr_std_matrices)
    Ah = matrix_index(A)
    matrix_trace = zero
    iprim = 0
    do np = 1, bundle%groups_on_node  ! Loop over partitions on processor
       if (bundle%nm_nodgroup(np) > 0) then
          do i = 1, bundle%nm_nodgroup(np) ! Loop over atoms in partition
             iprim = iprim + 1
             sf1 = mat(np,Ah)%ndimi(i)
             loc = (sf1 - 1) * sf1 + sf1 - 1 + mat(np,Ah)%onsite(i)
             if (loc > mat_p(A)%length) &
                  call cq_abort("Overflow error in trace: ", &
                                loc, mat_p(A)%length)
             do j = 1, sf1
                loc = (j - 1) * sf1 + j - 1 + mat(np,Ah)%onsite(i)
                matrix_trace = matrix_trace + mat_p(A)%matrix(loc) ! Add
             end do  ! End loop over NSF
          end do ! End loop over atoms in partition
       end if ! End if atoms in partition
    end do ! End loop over partitions on processor
    call gsum(matrix_trace) ! Global sum
    ! call stop_timer(tmr_std_matrices)
  end function matrix_trace
  !!***


  !! Atom-by-atom matrix_product_trace e.g. for atom charges
  !!   Call as atom_trace(A,B,prim%n_prim,res) where A and B are
  !!   matrix labels (e.g. matK, matS), len is bundle%n_prim and res
  !!   is a real(double) array of size len passed in.
  !! D.R. Bowler, 2008/03/21
  !!
  subroutine atom_trace(A, B, len, res)

    use datatypes
    use numbers,        only: zero
    use matrix_data,    only: mat
    use primary_module, only: bundle
    use GenBlas,        only: dot

    implicit none

    ! Input variables
    integer :: A,B, len
    real(double), dimension(len) :: res

    ! Local variables
    integer :: iprim, np, i, l, Ah, ist

!    call start_timer(tmr_std_matrices)
    res = zero
    Ah = matrix_index(A)
    iprim = 0
    if (matrix_index(A) == matrix_index(B)) then ! The matrices are the same range
       do np = 1, bundle%groups_on_node ! Loop over partitions on processor
          if (bundle%nm_nodgroup(np) > 0) then ! if there are atoms in partition
             do i = 1, bundle%nm_nodgroup(np) ! loop over atoms in partition
                iprim = iprim + 1   ! count no. of atoms
                ist = mat(np,Ah)%nd_offset + mat(np,Ah)%i_nd_acc(i)
                if (i < bundle%nm_nodgroup(np)) then 
                   ! if still going through partition atoms:
                   l = mat(np,Ah)%i_nd_acc(i+1) - mat(np,Ah)%i_nd_acc(i)
                else ! if finished all atoms in partition
                   l = mat(np,Ah)%part_nd_nabs - mat(np,Ah)%i_nd_acc(i) + 1
                end if
                ! atomic charge for iprim atom:
                res(iprim) = dot(l, mat_p(A)%matrix(ist:ist+l-1), 1, &
                                 mat_p(B)%matrix(ist:ist+l-1),1)
             end do
          end if
       end do ! recovered the charges per processor, need to gsum
    end if
    ! call stop_timer(tmr_std_matrices)

    return
  end subroutine atom_trace
  !!***


  !!****f* mult_module/matrix_product_trace
  !! PURPOSE
  !!   Calculate Tr(AB)
  !! INPUTS
  !! OUTPUT
  !! RETURN VALUE
  !! AUTHOR
  !!   D. Bowler
  !! CREATION DATE 
  !!   
  !! MODIFICATION HISTORY
  !!   2012/02/08 L.Tong
  !!   - Added intent(in) attribute for A and B.
  !!   2014/09/22 lat
  !!   - Added printing of matrix length  
  !! SOURCE
  !!
  function matrix_product_trace(A, B)

    use datatypes
    use numbers
    use GenBlas,     only: dot, vdot, asum

    implicit none

    ! results
    real(double) :: matrix_product_trace
    
    ! Passed variables
    integer, intent(in) :: A, B

    ! call start_timer(tmr_std_matrices)
    ! If the matrices are the same range
    if (matrix_index(A) == matrix_index(B)) then
       if (mat_p(A)%length /= mat_p(B)%length) &
            call cq_abort("Length error in mat_prod_tr: ", &
                          mat_p(A)%length, mat_p(B)%length)
       matrix_product_trace = &
            vdot(mat_p(A)%length, mat_p(A)%matrix, 1, mat_p(B)%matrix, 1)
    else
       ! For now I'm going to abort; we could have (say) optional
       ! argument giving the mult-type and range of the product
       call cq_abort("Matrix product for two different ranges not implemented", &
            mat_p(A)%length, mat_p(B)%length)
    end if
    ! call stop_timer(tmr_std_matrices)

    return
  end function matrix_product_trace
  !!*****


  !!****f* mult_module/matrix_product_trace_length
  !! PURPOSE
  !!   Calculate Tr(AB)
  !! INPUTS
  !! OUTPUT
  !! RETURN VALUE
  !! AUTHOR
  !!   D. Bowler/L. Truflandier
  !! CREATION DATE 
  !!   
  !! MODIFICATION HISTORY
  !!   2014/09/22 lat
  !!   - Added printing of matrix length  
  !!   - delete check with respect to matrix_index
  !! SOURCE
  !!
  function matrix_product_trace_length(A, B)

    use datatypes
    use numbers
    use GenBlas,     only: dot, vdot, asum

    implicit none

    ! results
    real(double) :: matrix_product_trace_length
    
    ! Passed variables
    integer, intent(in) :: A, B

    ! call start_timer(tmr_std_matrices)
    ! If the matrices are the same range
    !if (matrix_index(A) == matrix_index(B)) then
    if (mat_p(A)%length /= mat_p(B)%length) then
       call cq_abort("Length error in mat_prod_tr: ", &
            mat_p(A)%length, mat_p(B)%length)
    else
       matrix_product_trace_length = &
         vdot(mat_p(A)%length, mat_p(A)%matrix, 1, mat_p(B)%matrix, 1)
    end if
    !else
    ! For now I'm going to abort; we could have (say) optional
    ! argument giving the mult-type and range of the product
    !   call cq_abort("Matrix product for two different ranges not implemented", &
    !        mat_p(A)%length, mat_p(B)%length)

    ! call stop_timer(tmr_std_matrices)

    return
  end function matrix_product_trace_length
!!*****


  function return_matrix_len(A)
    implicit none
    ! Passed variable
    integer :: return_matrix_len
    integer :: A
    return_matrix_len = mat_p(A)%length
  end function return_matrix_len

  ! i is the primary set number
  ! j is gcs%icover_ibeg(np)+ni - 1
  function return_matrix_value(A, np, ni, ip, nabj, f1, f2, onsite)

    use datatypes
    use basic_types,   only: primary_set, cover_set
    use matrix_module, only: matrix_halo, matrix
    use matrix_data,   only: halo, mat
    use cover_module,  only: BCS_parts

    implicit none

    ! result
    real(double) :: return_matrix_value
    ! Passed variables
    integer :: A, np, ni, ip, nabj, f1, f2
    integer, optional :: onsite
    ! Local variables
    integer :: Ah, pos, jnp, jni, j, sf1

!    call start_timer(tmr_std_matrices)
    Ah = matrix_index(A)
    sf1 = mat(np,Ah)%ndimi(ni)
    if(present(onsite)) then
       pos = mat(np,Ah)%onsite(ni)
    else
       jnp = mat(np,Ah)%i_part(mat(np,Ah)%i_acc(ni)+nabj-1)
       jni = mat(np,Ah)%i_seq(mat(np,Ah)%i_acc(ni)+nabj-1)
       j = BCS_parts%icover_ibeg(jnp)+jni-1
       pos = halo(Ah)%i_h2d( (ip-1)*halo(Ah)%ni_in_halo + halo(Ah)%i_halo(j) )
    end if
    ! In future we'll just have matrix_pos here
    if ((f2-1)*sf1+f1-1+pos > mat_p(A)%length) &
         call cq_abort("Error in return_matrix_value", &
                       (f2-1)*sf1+f1-1+pos,mat_p(A)%length)
    return_matrix_value = mat_p(A)%matrix((f2-1)*sf1+f1-1+pos)
    ! call stop_timer(tmr_std_matrices)
    !if(pos/=posin) write(io_lun,*) 'Pos: ',pos,posin,i,j,halo(Ah)%i_halo(j),Ah
  end function return_matrix_value
!!*****


  function return_matrix_value_pos(A, i)
    use datatypes
    use GenComms, only: cq_abort

    implicit none

    ! result
    real(double) :: return_matrix_value_pos
    ! Passed variables
    integer, intent(in) :: A, i

!    call start_timer(tmr_std_matrices)
    if (i > mat_p(A)%length) &
         call cq_abort("Matrix overrun in return_matrix_value_pos: ", i, A)
    ! In future we'll just have matrix_pos here
    return_matrix_value_pos = mat_p(A)%matrix(i)
  end function return_matrix_value_pos
!!*****


  subroutine return_matrix_block_pos(A, i, val, size)
    use datatypes
    use GenComms, only: cq_abort

    implicit none

    ! Passed variables
    integer :: A, i,size
    real(double), dimension(size) :: val

    ! call start_timer(tmr_std_matrices)
    if (i+size-1 > mat_p(A)%length) &
         call cq_abort("Matrix overrun in return_matrix_value_pos: ", i, A)
    ! In future we'll just have matrix_pos here
    val(1:size) = mat_p(A)%matrix(i:i+size-1)
    ! call stop_timer(tmr_std_matrices)
  end subroutine return_matrix_block_pos
!!*****


  function matrix_pos(A, iprim, j_in_halo, f1, f2)

    use matrix_module, only: matrix_halo
    use matrix_data,   only: halo
    use GenComms,      only: cq_abort

    implicit none

    ! Result
    integer  :: matrix_pos
    ! Passed variables
    integer, intent(in) :: A, iprim, j_in_halo
    integer, optional, intent(in) :: f1, f2
    ! Local variables
    integer :: Ah, sf1

!    call start_timer(tmr_std_matrices)
    if (A < 1) call cq_abort("Error in matrix_pos: matrix number illegal: ", A)
    Ah = matrix_index(A)

    if (j_in_halo > halo(Ah)%mx_halo) &
         call cq_abort('Halo indexing error in matrix_pos: ', &
                       j_in_halo, halo(Ah)%mx_halo)
    if (present(f1) .and. present(f2)) then
       sf1 = halo(Ah)%ndimi(iprim)
       matrix_pos = (f2 - 1) * sf1 + f1 - 1 + &
                    halo(Ah)%i_h2d((iprim-1)*halo(Ah)%ni_in_halo + j_in_halo)
    else
       matrix_pos = halo(Ah)%i_h2d((iprim-1)*halo(Ah)%ni_in_halo + j_in_halo)
    end if
    ! Some length check here ?
    if (matrix_pos > mat_p(A)%length) &
         call cq_abort("Overrun error in matrix_pos: ", &
                       matrix_pos, mat_p(A)%length)
    ! call stop_timer(tmr_std_matrices)
  end function matrix_pos

  subroutine scale_matrix_value(A, np, ni, ip, nabj, f1, f2, val, &
                                onsite)
    use datatypes
    use matrix_module, only: matrix_halo, matrix
    use matrix_data,   only: halo, mat
    use cover_module,  only: BCS_parts

    implicit none

    ! Passed variables
    integer :: A, np, ni, ip, nabj, f1, f2
    integer, optional :: onsite
    real(double) :: val
    ! Local variables
    integer :: Ah, pos, j, jnp, jni, sf1

!    call start_timer(tmr_std_matrices)
    Ah = matrix_index(A)
    sf1 = halo(Ah)%ndimi(ip)
    if (present(onsite)) then
       pos = mat(np,Ah)%onsite(ni)
    else
       jnp = mat(np,Ah)%i_part(mat(np,Ah)%i_acc(ni)+nabj-1)
       jni = mat(np,Ah)%i_seq(mat(np,Ah)%i_acc(ni)+nabj-1)
       j = BCS_parts%icover_ibeg(jnp)+jni-1
       pos = halo(Ah)%i_h2d( (ip-1)*halo(Ah)%ni_in_halo + halo(Ah)%i_halo(j) )
    end if
    ! In future we'll just have matrix_pos here
    if (f1-1+(f2-1)*sf1+pos > mat_p(A)%length) &
         call cq_abort("Overrun error in scale_matrix_value: ", &
                       f1-1+(f2-1)*sf1+pos, mat_p(A)%length)
    mat_p(A)%matrix(f1-1+(f2-1)*sf1+pos) = mat_p(A)%matrix(f1-1+(f2-1)*sf1+pos)*val
!    call stop_timer(tmr_std_matrices)

    return
  end subroutine scale_matrix_value
!!*****


  subroutine store_matrix_value(A, np, ni, ip, nabj, f1, f2, val, &
                                onsite)
    use datatypes
    use matrix_module, only: matrix_halo, matrix
    use matrix_data,   only: halo, mat
    use cover_module,  only: BCS_parts
    use GenComms,      only: cq_abort

    implicit none

    ! Passed variables
    integer :: A, np, ni, ip, nabj, f1, f2
    integer, optional :: onsite
    real(double) :: val
    ! Local variables
    integer :: Ah, pos, j, jnp, jni, sf1

!    call start_timer(tmr_std_matrices)
    Ah = matrix_index(A)
    sf1 = halo(Ah)%ndimi(ip)
    if(present(onsite)) then
       pos = mat(np,Ah)%onsite(ni)
    else
       jnp = mat(np,Ah)%i_part(mat(np,Ah)%i_acc(ni)+nabj-1)
       jni = mat(np,Ah)%i_seq(mat(np,Ah)%i_acc(ni)+nabj-1)
       j = BCS_parts%icover_ibeg(jnp)+jni-1
       pos = halo(Ah)%i_h2d((ip-1)*halo(Ah)%ni_in_halo + halo(Ah)%i_halo(j))
    end if
    if(f1-1+(f2-1)*sf1+pos > mat_p(A)%length) &
         call cq_abort("Overrun error in store_matrix_value: ", &
                       f1-1+(f2-1)*sf1+pos, mat_p(A)%length)
    ! In future we'll just have matrix_pos here
    mat_p(A)%matrix(f1-1+(f2-1)*sf1+pos) = &
         mat_p(A)%matrix(f1-1+(f2-1)*sf1+pos) + val
!    call stop_timer(tmr_std_matrices)
  end subroutine store_matrix_value
!!*****


  subroutine store_matrix_value_pos(A, i, val)

    use datatypes
    use GenComms, only: cq_abort

    implicit none

    ! Passed variables
    integer,      intent(in) :: A, i
    real(double), intent(in) :: val

!    call start_timer(tmr_std_matrices)
    if (A < 1) &
         call cq_abort("Error in matrix_pos: matrix number illegal: ", A)
    if (i > mat_p(A)%length .or. i < 1) &
         call cq_abort("Matrix overrun in store_matrix_value_pos: ", &
                       i, mat_p(A)%length)
    ! if (abs(val)>100.0_double) write(io_lun,*) 'Error: ',A,i,val
    ! In future we'll just have matrix_pos here
    mat_p(A)%matrix(i) = mat_p(A)%matrix(i) + val
    ! call stop_timer(tmr_std_matrices)
    ! if (abs(mat_p(A)%matrix(i))>100.0_double) write(io_lun,*) 'Error: ',A,i,val
   
    return
  end subroutine store_matrix_value_pos
!!*****


  subroutine store_matrix_block_pos(A, i, val, size)

    use datatypes
    use GenComms, only: cq_abort
    use GenBlas,  only: axpy
    use numbers,  only: one

    implicit none

    ! Passed variables
    integer :: A, i, size
    real(double) :: val(size)

    ! call start_timer(tmr_std_matrices)
    if (i+size-1 > mat_p(A)%length) &
         call cq_abort("Matrix overrun in store_matrix_value_pos: ", &
                       i+size-1, mat_p(A)%length)
    call axpy(size, one, val, 1, mat_p(A)%matrix(i:), 1)
    ! call stop_timer(tmr_std_matrices)

    return
  end subroutine store_matrix_block_pos
!!*****


  ! subroutine store_onsite_matrix_value(A,np,i,f1,f2,val)
  !
  !   use datatypes
  !   use matrix_module, only: matrix
  !   use matrix_data, only: mat
  !
  !   implicit none
  !
  !   ! Passed variables
  !   integer :: A, np, i, f1, f2
  !   real(double) :: val
  !   ! Local variables
  !   integer :: Ah, pos
  !
  !   Ah = matrix_index(A)
  !   pos = mat(np,Ah)%onsite(i)
  !   ! In future we'll just have matrix_pos here
  !   mat_p(A)%matrix(f1,f2,pos) = val
  ! end subroutine store_onsite_matrix_value

  !!****f* mult_module/allocate_temp_matrix *
  !!
  !!  NAME
  !!   allocate_temp_matrix
  !!  USAGE
  !!  PURPOSE
  !!   Allocates space for temporary matrix and zeroes
  !!  INPUTS
  !!  USES
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   26/04/00
  !!  MODIFICATION HISTORY
  !!   2016/05/09 10:58 dave
  !!    Added output to help user increase number of temporary matrices
  !!  SOURCE
  !!
  function allocate_temp_matrix(range, trans_range, f1_type, f2_type)

    use matrix_data,   only: mat
    use GenComms,      only: cq_abort, myid
    use numbers,       only: zero
    use global_module, only: iprint_mat, area_matrices
    use GenBlas,       only: scal
    use memory_module, only: reg_alloc_mem, type_dbl

    implicit none
    
    ! Result
    integer :: allocate_temp_matrix
    ! Passed variables
    integer :: range, trans_range
    integer, optional :: f1_type, f2_type
    ! Local variables
    integer :: stat, i

    if (current_matrix == max_matrices) then
       if(myid==0) write(io_lun,fmt='(4x,"You have exceeded the maximum number of temporary matrices in Conquest.")') 
       if(myid==0) write(io_lun,fmt='(4x,"Consider increasing the limit on temporary matrices with General.MaxTempMatrices")') 
       call cq_abort("Overrun number of matrices ", max_matrices)
    end if
    current_matrix = current_matrix + 1
    ! Identify the matrix
    if (range > mx_matrices) &
         call cq_abort("Array error in alloc_temp_mat: ", range, mx_matrices)
    matrix_index(current_matrix) = range
    trans_index(current_matrix) = trans_range
    if (present(f1_type) .and. present(f2_type)) then
       mat_p(current_matrix)%sf1_type = f1_type
       mat_p(current_matrix)%sf2_type = f2_type
    else
       mat_p(current_matrix)%sf1_type = sf
       mat_p(current_matrix)%sf2_type = sf
    end if
    mat_p(current_matrix)%length = mat(1,range)%length
    call start_timer(tmr_std_allocation)
    allocate(mat_p(current_matrix)%matrix(mat_p(current_matrix)%length), &
             STAT=stat)
    if (stat /= 0) call cq_abort("Error allocating matrix ", &
                                 current_matrix, mat_p(current_matrix)%length)
    call stop_timer(tmr_std_allocation)
    ! call scal(mat_p(current_matrix)%length,zero,mat_p(current_matrix)%matrix,1)
    do i = 1, mat_p(current_matrix)%length
       mat_p(current_matrix)%matrix(i) = zero
    end do
    allocate_temp_matrix = current_matrix
    if (iprint_mat > 4 .and. myid == 0) &
         write(io_lun,fmt='(10x,a,2i5)') '(A)Current matrix is: ', &
                         current_matrix, mat_p(current_matrix)%length
    call reg_alloc_mem(area_matrices, mat_p(current_matrix)%length, type_dbl)

    return
  end function allocate_temp_matrix
!!***


  !!****f* mult_module/free_temp_matrix *
  !!
  !!  NAME
  !!   free_temp_matrix
  !!  USAGE
  !!  PURPOSE
  !!   Frees space for temporary matrix and zeroes
  !!  INPUTS
  !!  USES
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   26/04/00
  !!  MODIFICATION HISTORY
  !!    
  !!  SOURCE
  !!
  subroutine free_temp_matrix(A)

    use GenComms,      only: cq_abort, myid
    use global_module, only: iprint_mat, area_matrices
    use memory_module, only: reg_dealloc_mem, type_dbl

    implicit none

    ! Passed variables
    integer :: A
    ! Local variables
    integer :: stat

    call reg_dealloc_mem(area_matrices, mat_p(current_matrix)%length, &
                         type_dbl)
    if (A /= current_matrix) &
         call cq_abort("Out-of-order deallocation of matrices ", &
                       A, current_matrix)
    if (iprint_mat > 4 .and. myid == 0) &
         write (io_lun,fmt='(10x,a,2i5)') '(D)Current matrix is: ', &
                          current_matrix, mat_p(current_matrix)%length
    call start_timer(tmr_std_allocation)
    deallocate(mat_p(current_matrix)%matrix, STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error deallocating matrix ", &
                       current_matrix, mat_p(current_matrix)%length)
    call stop_timer(tmr_std_allocation)
    mat_p(current_matrix)%length = 0
    matrix_index(current_matrix) = 0
    current_matrix = current_matrix - 1

    return
  end subroutine free_temp_matrix
!!***

  !!****f* mult_module/AtomF_to_SF_transform *
  !!
  !!  NAME
  !!   AtomF_to_SF_transform
  !!  USAGE
  !!  PURPOSE
  !!   Multiply SF coefficients to transform matAtomF -> matSF
  !!  INPUTS
  !!  USES
  !!  AUTHOR
  !!   A.Nakata
  !!  CREATION DATE
  !!   26/09/16
  !!  MODIFICATION HISTORY
  !!
  !!  SOURCE
  !!
  subroutine AtomF_to_SF_transform(matSF, matAtomF, spin, Arange)

    use numbers
    use GenComms,      only: cq_abort
    use global_module, only: flag_SpinDependentSF
    use matrix_data,   only: Srange, Hrange, aSs_range, aHs_range

    implicit none

    ! Passed variables
    integer :: matSF, matAtomF, spin, Arange
    ! Local variables
    integer :: spin_SF, mat_AtomF_SF, ind1_mult, ind2_mult


    spin_SF = 1
    if (flag_SpinDependentSF) spin_SF = spin

    if (Arange.eq.Srange) then
       mat_AtomF_SF = allocate_temp_matrix(aSs_range,0,atomf,sf)
       ind1_mult = aSa_sCaTr_aSs
       ind2_mult = sCa_aSs_sSs
    else if (Arange.eq.Hrange) then
       mat_AtomF_SF = allocate_temp_matrix(aHs_range,0,atomf,sf)
       ind1_mult = aHa_sCaTr_aHs
       ind2_mult = sCa_aHs_sHs
    else
       call cq_abort("AtomF_to_SF_transform supports only Srange and Hrange.")
    endif

    call matrix_scale(zero, mat_AtomF_SF)
    call matrix_scale(zero, matSF)

    ! Transform matAatomf(atomf,atomf) in Aatomf_range -> matA(sf,sf) in Arange
    ! step 1: B(atomf,sf) = A(atomf,atomf) * SFcoeff(sf,atomf)-DAGGER
    call matrix_product(matAtomF, matSFcoeff_tran(spin_SF), mat_AtomF_SF, mult(ind1_mult))
    ! step 2: A(sf,sf) = SFcoeff(sf,atomf) * B(atomf,sf)
    call matrix_product(matSFcoeff(spin_SF), mat_AtomF_SF, matSF, mult(ind2_mult) )

    call free_temp_matrix(mat_AtomF_SF)
!
    return
  end subroutine AtomF_to_SF_transform
!!***

  !!****f* mult_module/SF_to_AtomF_transform *
  !!
  !!  NAME
  !!   SF_to_AtomF_transform
  !!  USAGE
  !!  PURPOSE
  !!   Multiply SF coefficients to transform matSF -> matAtomF
  !!  INPUTS
  !!  USES
  !!  AUTHOR
  !!   A.Nakata
  !!  CREATION DATE
  !!   26/09/16
  !!  MODIFICATION HISTORY
  !!
  !!  SOURCE
  !!
  subroutine SF_to_AtomF_transform(matSF, matAtomF, spin, Arange)

    use numbers
    use GenComms,      only: cq_abort
    use global_module, only: flag_SpinDependentSF
    use matrix_data,   only: Srange, Hrange, aSs_range, aHs_range

    implicit none

    ! Passed variables
    integer :: matSF, matAtomF, spin, Arange
    ! Local variables
    integer :: spin_SF, mat_AtomF_SF, ind1_mult, ind2_mult


    spin_SF = 1
    if (flag_SpinDependentSF) spin_SF = spin

    if (Arange.eq.Srange) then
       mat_AtomF_SF = allocate_temp_matrix(aSs_range,0,atomf,sf)
       ind1_mult = sCaTr_sSs_aSs
       ind2_mult = aSs_sCa_aSa
    else if (Arange.eq.Hrange) then
       mat_AtomF_SF = allocate_temp_matrix(aHs_range,0,atomf,sf)
       ind1_mult = sCaTr_sHs_aHs
       ind2_mult = aHs_sCa_aHa
    else
       call cq_abort("SF_to_AtomF_transform supports only Srange and Hrange.")
    endif

    call matrix_scale(zero, mat_AtomF_SF)
    call matrix_scale(zero, matAtomF)

    ! Transform matA(sf,sf) in Arange -> matAatomf(atomf,atomf) in Aatomf_range
    ! step 1: B(atomf,sf) = SFcoeff(sf,atomf)-DAGGER * A(sf,sf)
    call matrix_product(matSFcoeff_tran(spin_SF), matSF, mat_AtomF_SF, mult(ind1_mult))
    ! step 2: A(atomf,atomf) = B(atomf,sf) * SFcoeff(sf,atomf) 
    !         use SFcoeff_tran instead of SFcoeff because mult_type=2
    call matrix_product(mat_AtomF_SF, matSFcoeff_tran(spin_SF), matAtomF, mult(ind2_mult))

    call free_temp_matrix(mat_AtomF_SF)
!
    return
  end subroutine SF_to_AtomF_transform
!!***

!!****f* multisiteSF_module/check_InterAtomicDistances *
!!
!!  NAME
!!   check_InterAtomicDistances
!!  USAGE
!!
!!  PURPOSE
!!   Check interatomic distances in given coordinates
!!   whetehr the atoms have appropiate distances (not too close to each other),
!!   using the information of matSatomf.
!!   Ghost atoms will be ignored.
!!
!!  AUTHOR
!!   A.Nakata
!!  CREATION DATE
!!   2019/02/22
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine check_InterAtomicDistances

    use datatypes
    use global_module,  ONLY: id_glob, species_glob
    use species_module, ONLY: type_species
    use GenComms,       ONLY: cq_abort
    use group_module,   ONLY: parts
    use primary_module, ONLY: bundle
    use cover_module,   ONLY: BCS_parts
    use matrix_data,    ONLY: mat, Srange

    implicit none

    ! Local variables
    integer :: iprim, np, i, j, ist, atom_num, atom_i, atom_spec
    integer :: neigh_global_part, atom_j, neigh_species
    integer :: gcspart
    real(double) :: dx, dy, dz, r2, r_min

    r_min = 100.0_double

    iprim = 0
    do np = 1,bundle%groups_on_node ! Loop over primary set partitions
       if(bundle%nm_nodgroup(np)>0) then ! If there are atoms in partition
          do i = 1,bundle%nm_nodgroup(np) ! Loop over atoms
             atom_num = bundle%nm_nodbeg(np) + i - 1
             iprim = iprim + 1
             atom_i = bundle%ig_prim(iprim)   ! global number of i
             atom_spec = bundle%species(atom_num)
             if(type_species(atom_spec) < 0) cycle   ! if ghost atoms
             do j = 1, mat(np,Srange)%n_nab(i) ! Loop over neighbours of atom
                ist = mat(np,Srange)%i_acc(i) + j - 1
                ! Build the distances between atoms
                gcspart = BCS_parts%icover_ibeg(mat(np,Srange)%i_part(ist))+mat(np,Srange)%i_seq(ist)-1
                neigh_global_part = BCS_parts%lab_cell(mat(np,Srange)%i_part(ist))
                ! global number of j
                atom_j = id_glob(parts%icell_beg(neigh_global_part)+mat(np,Srange)%i_seq(ist)-1) 
                neigh_species = species_glob(atom_j)
                if(type_species(neigh_species) < 0) cycle   ! if ghost atoms
                if (atom_i/=atom_j .AND. mat(np,Srange)%radius(ist)<=InitAtomicDistance_min) then
                   call cq_abort("Atoms are too close to each other:",atom_i,atom_j)
                endif
                if (mat(np,Srange)%radius(ist)<r_min) r_min = mat(np,Srange)%radius(ist)
             end do ! j
          end do ! i
       end if ! End if nm_nodgroup > 0
    end do ! np
!
    if (r_min>InitAtomicDistance_max) call cq_abort("Atoms are too far apart! Minimum distance is ",r_min)
!
    return
  end subroutine check_InterAtomicDistances
  !!***

end module mult_module
