! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
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
!!    Made all long-range (LH or more) matrices allocatable in main_matrix_multiply and 
!!    improved headers slightly.  Also added RCS static object.
!!   12:29, 04/02/2003 drb 
!!    Added Trace
!!   13:51, 10/02/2003 drb 
!!    Added end_ops and fmmi to deallocate at the end of a run
!!   14:45, 26/02/2003 drb 
!!    Added implicit none
!!   10:42, 06/03/2003 drb 
!!    Changed prim to bundle in main_matrix_multiply (ifc problem with aliasing)
!!   14:28, 16/11/2004 dave 
!!    Removed TS_TS_T multiply (problematic !)
!!   10:58, 13/02/2006 drb 
!!    Complete new style of matrix manipulation which hides matrix data from rest of code: new routines, 
!!    variables etc.
!!   2008/02/06 08:25 dave
!!    Changed for output to file not stdout
!!  TODO
!!   17/06/2002 dave
!!    Split main_matrix_multiply into sensible small routines (like getK, getDOmega)
!!    Change immi so that the multiplication initialisation is done by a subroutine and call it many times
!!   2008/05/22 ast
!!    Added timers
!!  SOURCE
!!
module mult_module

  use datatypes
  use matrix_module
  use matrix_data, ONLY: mx_matrices, matrix_pointer
  use global_module, ONLY: sf, nlpf, paof, io_lun
  use GenComms, ONLY: cq_abort
  use timer_stdclocks_module, ONLY: start_timer,stop_timer,tmr_std_allocation

  implicit none
  save

  integer(integ), parameter :: H_SP_SP  = 1   ! type 2
  integer(integ), parameter :: L_S_LS   = 2   ! type 1
  integer(integ), parameter :: L_H_LH   = 3   ! type 1
  integer(integ), parameter :: T_S_TS   = 4   ! type 1
!  integer(integ), parameter :: S_TS_S   = 5   ! type 2, maybe do ST_S_S
!  integer(integ), parameter :: S_LS_L   = 6   ! type 2, maybe SL_S_L (prune 7?)
  integer(integ), parameter :: SL_S_SLS = 7   ! type 1
  integer(integ), parameter :: LH_L_SLS = 8   ! type 1 or 2 depending on H
  integer(integ), parameter :: LS_L_LSL = 9   ! type 1
  integer(integ), parameter :: HL_S_LSL = 10  ! type 1 or 2 depending on H
  integer(integ), parameter :: HLS_LS_L = 11  ! type 2
  integer(integ), parameter :: LSL_SL_L = 11  ! type 2
  integer(integ), parameter :: SLS_LS_L = 12  ! type 2
  integer(integ), parameter :: LHL_SL_S = 13  ! type 2
  integer(integ), parameter :: LSL_SL_H = 14  ! type 2
  integer(integ), parameter :: SP_PS_H  = 15  ! type 1 (possibly 2)
  integer(integ), parameter :: TS_T_T   = 16  ! type 2
  integer(integ), parameter :: TH_T_L   = 17  ! type 2 or 1 depending on H
  integer(integ), parameter :: T_H_TH   = 18  ! type 1
  integer(integ), parameter :: T_L_TL   = 19  ! type 1
  integer(integ), parameter :: TL_T_L   = 20  ! type 2
  integer(integ), parameter :: PAOP_PS_H   = 21  ! type 1/2 (cf SP_PS_H)

  integer(integ), parameter :: mx_mults = 21

  type(matrix_mult) :: mult(mx_mults)

  integer(integ), parameter :: S_trans   = 1
  integer(integ), parameter :: L_trans   = 2
  integer(integ), parameter :: T_trans   = 3
  integer(integ), parameter :: SP_trans  = 4 
  integer(integ), parameter :: LS_trans  = 5
  integer(integ), parameter :: LH_trans  = 6
  integer(integ), parameter :: LSL_trans = 7

  integer(integ), parameter :: mx_trans = 7

  type(pair_data), allocatable, dimension(:,:) :: pairs
  integer, dimension(:), pointer :: Spairind, Lpairind, Tpairind, SPpairind, LSpairind, LHpairind, LSLpairind

  type(matrix_trans), target :: ltrans(mx_matrices)
  type(trans_remote) :: gtrans(mx_trans)

  type(matrix_pointer), allocatable, dimension(:) :: mat_p
  integer :: max_matrices, current_matrix
  integer, public :: matS, matH, matL, matT, matTtran, matLS, matSL, matSC, matCS, matK, matM4, matU, matphi, matM12, &
       matUT, matKE, matNL, matdH, matdS
  integer, allocatable, dimension(:), public :: matrix_index, trans_index

  ! -------------------------------------------------------
  ! RCS ident string for object file id
  ! -------------------------------------------------------
  character(len=80), private :: RCSid = "$Id$"
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
!!  SOURCE
!!
  subroutine immi(parts,prim,gcs,myid,partial)

    use datatypes
    use basic_types
    use matrix_data
    use matrix_elements_module
    use mult_init_module
    use maxima_module, ONLY: maxpartsproc, maxnabaprocs
    use GenComms, ONLY: my_barrier, cq_abort
    use matrix_comms_module, ONLY: find_neighbour_procs
    use global_module, ONLY: area_matrices
    use memory_module, ONLY: reg_alloc_mem, type_int

    implicit none

    ! Passed variables
    type(group_set), target :: parts
    type(primary_set), target :: prim
    type(cover_set), target :: gcs
    integer :: myid
    integer, OPTIONAL :: partial

    ! Local variables
    integer :: stat,i,j
    real(double) :: ra,rc

!    call start_timer(tmr_std_matrices)
    max_matrices = mx_matrices + 100 ! Leave plenty of room for temporary matrices
    if(.not.allocated(mat_p)) then
       call start_timer(tmr_std_allocation)
       allocate(mat_p(max_matrices), matrix_index(max_matrices), trans_index(max_matrices), STAT=stat)
       if(stat/=0) call cq_abort("Error allocating matrix indexing ",max_matrices,stat)
       call reg_alloc_mem(area_matrices,2*max_matrices,type_int)
       call stop_timer(tmr_std_allocation)
       do i=1,max_matrices
          mat_p(i)%length = 0
          nullify(mat_p(i)%matrix)
          matrix_index(i) = 0
          trans_index(i) = 0
       end do
    end if    
    if(.not.allocated(mat)) then
       call start_timer(tmr_std_allocation)
       allocate(mat(maxpartsproc,mx_matrices), STAT=stat)
       if(stat/=0) call cq_abort("Error allocating matrix type ",mx_matrices,maxpartsproc)
       call stop_timer(tmr_std_allocation)
    end if
    if(.not.allocated(halo)) then
       call start_timer(tmr_std_allocation)
       allocate(halo(mx_matrices), STAT=stat)
       if(stat/=0) call cq_abort("Error allocating halo type ",mx_matrices,stat)
       call stop_timer(tmr_std_allocation)
    end if
    mat%sf1_type = sf
    mat%sf2_type = sf
    ! First, initialise matrix indexing for mults and transposes
    ! for the different ranges
    call matrix_ini(parts,prim,gcs,mat(1:prim%groups_on_node,Lrange),Lmatind, &
         rcut(Lrange),myid-1,halo(Lrange),ltrans(Lrange))
    call matrix_ini(parts,prim,gcs,mat(1:prim%groups_on_node,Srange),Smatind, &
         rcut(Srange),myid-1,halo(Srange),ltrans(Srange))
    call matrix_ini(parts,prim,gcs,mat(1:prim%groups_on_node,Hrange),Hmatind, &
         rcut(Hrange),myid-1,halo(Hrange),ltrans(Hrange))
    mat(1:prim%groups_on_node,SPrange)%sf1_type = sf
    mat(1:prim%groups_on_node,SPrange)%sf2_type = nlpf
    call matrix_ini(parts,prim,gcs,mat(1:prim%groups_on_node,SPrange), SPmatind, &
         rcut(SPrange),myid-1,halo(SPrange),ltrans(SPrange))
    mat(1:prim%groups_on_node,PSrange)%sf1_type = nlpf
    mat(1:prim%groups_on_node,PSrange)%sf2_type = sf
    call matrix_ini(parts,prim,gcs,mat(1:prim%groups_on_node,PSrange), &
         PSmatind, rcut(SPrange),myid-1,halo(PSrange),ltrans(PSrange))
    call matrix_ini(parts,prim,gcs,mat(1:prim%groups_on_node,LSrange), LSmatind, &
         rcut(LSrange),myid-1,halo(LSrange),ltrans(LSrange))
    call matrix_ini(parts,prim,gcs,mat(1:prim%groups_on_node,LHrange), LHmatind, &
         rcut(LHrange),myid-1,halo(LHrange),ltrans(LHrange))
    call matrix_ini(parts,prim,gcs,mat(1:prim%groups_on_node,HLrange), &
         HLmatind, &
         rcut(LHrange),myid-1,halo(HLrange),ltrans(HLrange))
    call matrix_ini(parts,prim,gcs,mat(1:prim%groups_on_node,LSLrange), LSLmatind, &
         rcut(LSLrange),myid-1,halo(LSLrange),ltrans(LSLrange))
    call matrix_ini(parts,prim,gcs,mat(1:prim%groups_on_node,SLSrange), SLSmatind, &
         rcut(SLSrange),myid-1,halo(SLSrange),ltrans(SLSrange))
    call matrix_ini(parts,prim,gcs,mat(1:prim%groups_on_node,Trange),Tmatind, &
         rcut(Trange),myid-1,halo(Trange),ltrans(Trange))
    ! Add global transpose for the type 2 mults T?_T_L ?
    call matrix_ini(parts,prim,gcs,mat(1:prim%groups_on_node,TSrange), TSmatind, &
         rcut(TSrange),myid-1,halo(TSrange),ltrans(TSrange))
    call matrix_ini(parts,prim,gcs,mat(1:prim%groups_on_node,THrange), THmatind, &
         rcut(THrange),myid-1,halo(THrange),ltrans(THrange))
    call matrix_ini(parts,prim,gcs,mat(1:prim%groups_on_node,TLrange), TLmatind, &
         rcut(TLrange),myid-1,halo(TLrange),ltrans(TLrange))
    ! Transposes
    call matrix_ini(parts,prim,gcs,mat(1:prim%groups_on_node,LTrrange),LTrmatind, &
         rcut(Lrange),myid-1,halo(LTrrange),ltrans(LTrrange))
    call matrix_ini(parts,prim,gcs,mat(1:prim%groups_on_node,SLrange), &
         SLmatind, rcut(LSrange),myid-1,halo(SLrange),ltrans(SLrange))
    call matrix_ini(parts,prim,gcs,mat(1:prim%groups_on_node,TTrrange),TTrmatind, &
         rcut(Trange),myid-1,halo(TTrrange),ltrans(TTrrange))
    ! Differentials
    mat(1:prim%groups_on_node,dSrange)%sf1_type = paof
    mat(1:prim%groups_on_node,dSrange)%sf2_type = sf
    call matrix_ini(parts,prim,gcs,mat(1:prim%groups_on_node,dSrange),dSmatind, &
         rcut(Srange),myid-1,halo(dSrange),ltrans(dSrange))
    mat(1:prim%groups_on_node,dHrange)%sf1_type = paof
    mat(1:prim%groups_on_node,dHrange)%sf2_type = sf
    call matrix_ini(parts,prim,gcs,mat(1:prim%groups_on_node,dHrange),dHmatind, &
         rcut(Hrange),myid-1,halo(dHrange),ltrans(dHrange))
    mat(1:prim%groups_on_node,PAOPrange)%sf1_type = paof
    mat(1:prim%groups_on_node,PAOPrange)%sf2_type = nlpf
    call matrix_ini(parts,prim,gcs,mat(1:prim%groups_on_node,PAOPrange), &
         PAOPmatind, rcut(SPrange),myid-1,halo(PAOPrange),ltrans(PAOPrange))
    !if(.NOT.PRESENT(partial)) then
       call associate_matrices 
    !endif
    call find_neighbour_procs(parts,halo(LSLrange))
    call start_timer(tmr_std_allocation)
    allocate(pairs(maxnabaprocs+1,mx_trans),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating pairs in immi: ",maxnabaprocs+1,mx_trans)
    call stop_timer(tmr_std_allocation)
    do i=1,mx_trans
       do j=1,maxnabaprocs+1
          nullify(pairs(j,i)%submat)
       end do
    end do
    call trans_ini(parts,prim,gcs,mat(1:prim%groups_on_node,Lrange),myid-1,halo(Lrange),halo(Lrange),&
         ltrans(Lrange),gtrans(L_trans),pairs(:,L_trans),Lpairind)
    call trans_ini(parts,prim,gcs,mat(1:prim%groups_on_node,Srange),myid-1,halo(Srange),halo(Srange),&
         ltrans(Srange),gtrans(S_trans),pairs(:,S_trans),Spairind)
    call trans_ini(parts,prim,gcs,mat(1:prim%groups_on_node,SPrange),myid-1,halo(SPrange),halo(PSrange),&
         ltrans(SPrange),gtrans(SP_trans),pairs(:,SP_trans),SPpairind)
    call trans_ini(parts,prim,gcs,mat(1:prim%groups_on_node,LSrange),myid-1,halo(LSrange),halo(LSrange),&
         ltrans(LSrange),gtrans(LS_trans),pairs(:,LS_trans),LSpairind)
    call trans_ini(parts,prim,gcs,mat(1:prim%groups_on_node,LHrange),myid-1,halo(LHrange),halo(LHrange),&
         ltrans(LHrange),gtrans(LH_trans),pairs(:,LH_trans),LHpairind)
    call trans_ini(parts,prim,gcs,mat(1:prim%groups_on_node,LSLrange),myid-1,halo(LSLrange),halo(LSLrange),&
         ltrans(LSLrange),gtrans(LSL_trans),pairs(:,LSL_trans),LSLpairind)
    call trans_ini(parts,prim,gcs,mat(1:prim%groups_on_node,Trange),myid-1,halo(Trange),halo(Trange),&
         ltrans(Trange),gtrans(T_trans),pairs(:,T_trans),Tpairind)

    ! Now initialise the matrix multiplications
    ! Somewhere, we need to set mult_type 1/2 - probably here !
    ! 1-5
    ! (a) We need to arrange a global transpose before a type 2 mult
    ! (b) mult_wrap is passed A,B,C  - and send mat_mult C,B,A for a 
    !     type 2 mult
    ! (c) HotInvS_mm does TS_T_S as type two, and sends mult_wrap 
    !     TS, T, S
    mult(H_SP_SP)%mult_type = 2
    mult(H_SP_SP)%amat => mat(1:prim%groups_on_node,SPrange)
    mult(H_SP_SP)%bmat => mat(1:prim%groups_on_node,PSrange)
    mult(H_SP_SP)%cmat => mat(1:prim%groups_on_node,Hrange)
    mult(H_SP_SP)%ahalo => halo(SPrange)
    mult(H_SP_SP)%chalo => halo(Hrange)
    mult(H_SP_SP)%ltrans => ltrans(SPrange)
    !mult(H_SP_SP)%bindex => PSmatind
    mult(H_SP_SP)%parts => parts
    mult(H_SP_SP)%prim => prim
    mult(H_SP_SP)%gcs => gcs
    call mult_ini(mult(H_SP_SP),PSmatind,myid-1,prim%n_prim,parts)
    mult(L_S_LS)%mult_type = 1
    mult(L_S_LS)%amat => mat(1:prim%groups_on_node,Lrange)
    mult(L_S_LS)%bmat => mat(1:prim%groups_on_node,Srange)
    mult(L_S_LS)%cmat => mat(1:prim%groups_on_node,LSrange)
    mult(L_S_LS)%ahalo => halo(Lrange)
    mult(L_S_LS)%chalo => halo(LSrange)
    mult(L_S_LS)%ltrans => ltrans(Lrange)
    !mult(L_S_LS)%bindex => Smatind
    mult(L_S_LS)%parts => parts
    mult(L_S_LS)%prim => prim
    mult(L_S_LS)%gcs => gcs
    call mult_ini(mult(L_S_LS),Smatind,myid-1,prim%n_prim,parts)
    mult(L_H_LH)%mult_type = 1
    mult(L_H_LH)%amat => mat(1:prim%groups_on_node,Lrange)
    mult(L_H_LH)%bmat => mat(1:prim%groups_on_node,Hrange)
    mult(L_H_LH)%cmat => mat(1:prim%groups_on_node,LHrange)
    mult(L_H_LH)%ahalo => halo(Lrange)
    mult(L_H_LH)%chalo => halo(LHrange)
    mult(L_H_LH)%ltrans => ltrans(Lrange)
    !mult(L_H_LH)%bindex => Hmatind
    mult(L_H_LH)%parts => parts
    mult(L_H_LH)%prim => prim
    mult(L_H_LH)%gcs => gcs
    call mult_ini(mult(L_H_LH),Hmatind,myid-1,prim%n_prim,parts)
    mult(T_S_TS)%mult_type = 1
    mult(T_S_TS)%amat => mat(1:prim%groups_on_node,Trange)
    mult(T_S_TS)%bmat => mat(1:prim%groups_on_node,Srange)
    mult(T_S_TS)%cmat => mat(1:prim%groups_on_node,TSrange)
    mult(T_S_TS)%ahalo => halo(Trange)
    mult(T_S_TS)%chalo => halo(TSrange)
    mult(T_S_TS)%ltrans => ltrans(Trange)
    !mult(T_S_TS)%bindex => Smatind
    mult(T_S_TS)%parts => parts
    mult(T_S_TS)%prim => prim
    mult(T_S_TS)%gcs => gcs
    call mult_ini(mult(T_S_TS),Smatind,myid-1,prim%n_prim,parts)
    ! Removed 2006/03/28 DRB
    ! Not used: if re-instated, will require ST index
    !mult(S_TS_S)%mult_type = 2
    !mult(S_TS_S)%amat => mat(1:prim%groups_on_node,Srange)
    !mult(S_TS_S)%bmat => mat(1:prim%groups_on_node,TSrange)
    !mult(S_TS_S)%cmat => mat(1:prim%groups_on_node,Srange)
    !mult(S_TS_S)%ahalo => halo(Srange)
    !mult(S_TS_S)%chalo => halo(Srange)
    !mult(S_TS_S)%ltrans => ltrans(Srange)
    !mult(S_TS_S)%bindex => TSmatind
    !mult(S_TS_S)%parts => parts
    !mult(S_TS_S)%prim => prim
    !mult(S_TS_S)%gcs => gcs
    !call mult_ini(mult(S_TS_S),myid-1,prim%n_prim,parts)
    ! 6-10
    !mult(S_LS_L)%mult_type = 2
    !mult(S_LS_L)%amat => mat(1:prim%groups_on_node,Lrange)
    !mult(S_LS_L)%bmat => mat(1:prim%groups_on_node,LSrange)
    !mult(S_LS_L)%cmat => mat(1:prim%groups_on_node,Srange)
    !mult(S_LS_L)%ahalo => halo(Lrange)
    !mult(S_LS_L)%chalo => halo(Srange)
    !mult(S_LS_L)%ltrans => ltrans(Lrange)
    !mult(S_LS_L)%bindex => LSmatind
    !mult(S_LS_L)%parts => parts
    !mult(S_LS_L)%prim => prim
    !mult(S_LS_L)%gcs => gcs
    !call mult_ini(mult(S_LS_L),myid-1,prim%n_prim,parts)
    mult(SL_S_SLS)%mult_type = 1
    mult(SL_S_SLS)%amat => mat(1:prim%groups_on_node,SLrange)
    mult(SL_S_SLS)%bmat => mat(1:prim%groups_on_node,Srange)
    mult(SL_S_SLS)%cmat => mat(1:prim%groups_on_node,SLSrange)
    mult(SL_S_SLS)%ahalo => halo(SLrange)
    mult(SL_S_SLS)%chalo => halo(SLSrange)
    mult(SL_S_SLS)%ltrans => ltrans(SLrange)
    !mult(SL_S_SLS)%bindex => Smatind
    mult(SL_S_SLS)%parts => parts
    mult(SL_S_SLS)%prim => prim
    mult(SL_S_SLS)%gcs => gcs
    call mult_ini(mult(SL_S_SLS),Smatind,myid-1,prim%n_prim,parts)
    rc = (2.0_double*rcut(Srange)+rcut(Lrange))
    ra = rcut(Lrange)+rcut(Hrange)
    if(rc>=ra) then
       mult(LH_L_SLS)%mult_type = 1
       mult(LH_L_SLS)%amat => mat(1:prim%groups_on_node,LHrange)
       mult(LH_L_SLS)%bmat => mat(1:prim%groups_on_node,Lrange)
       mult(LH_L_SLS)%cmat => mat(1:prim%groups_on_node,SLSrange)
       mult(LH_L_SLS)%ahalo => halo(LHrange)
       mult(LH_L_SLS)%chalo => halo(SLSrange)
       mult(LH_L_SLS)%ltrans => ltrans(LHrange)
       !mult(LH_L_SLS)%bindex => Lmatind
       mult(LH_L_SLS)%parts => parts
       mult(LH_L_SLS)%prim => prim
       mult(LH_L_SLS)%gcs => gcs
       call mult_ini(mult(LH_L_SLS),Lmatind,myid-1,prim%n_prim,parts)
    else
       mult(LH_L_SLS)%mult_type = 2
       mult(LH_L_SLS)%amat => mat(1:prim%groups_on_node,SLSrange)
       mult(LH_L_SLS)%bmat => mat(1:prim%groups_on_node,LTrrange)
       mult(LH_L_SLS)%cmat => mat(1:prim%groups_on_node,LHrange)
       mult(LH_L_SLS)%ahalo => halo(SLSrange)
       mult(LH_L_SLS)%chalo => halo(LHrange)
       mult(LH_L_SLS)%ltrans => ltrans(SLSrange)
       !mult(LH_L_SLS)%bindex => LTrmatind
       mult(LH_L_SLS)%parts => parts
       mult(LH_L_SLS)%prim => prim
       mult(LH_L_SLS)%gcs => gcs
       call mult_ini(mult(LH_L_SLS),LTrmatind,myid-1,prim%n_prim,parts)
    endif
    mult(LS_L_LSL)%mult_type = 1
    mult(LS_L_LSL)%amat => mat(1:prim%groups_on_node,LSrange)
    mult(LS_L_LSL)%bmat => mat(1:prim%groups_on_node,Lrange)
    mult(LS_L_LSL)%cmat => mat(1:prim%groups_on_node,LSLrange)
    mult(LS_L_LSL)%ahalo => halo(LSrange)
    mult(LS_L_LSL)%chalo => halo(LSLrange)
    mult(LS_L_LSL)%ltrans => ltrans(LSrange)
    !mult(LS_L_LSL)%bindex => Lmatind
    mult(LS_L_LSL)%parts => parts
    mult(LS_L_LSL)%prim => prim
    mult(LS_L_LSL)%gcs => gcs
    call mult_ini(mult(LS_L_LSL),Lmatind,myid-1,prim%n_prim,parts)
    !ra = (rcut(Hrange)+rcut(Lrange))
    !rc = 2.0_double*rcut(Lrange)+rcut(Srange)
    !if(rc>=ra) then
       mult(HL_S_LSL)%mult_type = 1
       mult(HL_S_LSL)%amat => mat(1:prim%groups_on_node,HLrange)
       mult(HL_S_LSL)%bmat => mat(1:prim%groups_on_node,Srange)
       mult(HL_S_LSL)%cmat => mat(1:prim%groups_on_node,LSLrange)
       mult(HL_S_LSL)%ahalo => halo(HLrange)
       mult(HL_S_LSL)%chalo => halo(LSLrange)
       mult(HL_S_LSL)%ltrans => ltrans(HLrange)
       !mult(HL_S_LSL)%bindex => Smatind
       mult(HL_S_LSL)%parts => parts
       mult(HL_S_LSL)%prim => prim
       mult(HL_S_LSL)%gcs => gcs
       call mult_ini(mult(HL_S_LSL),Smatind,myid-1,prim%n_prim,parts)
    !else
    !   mult(HL_S_LSL)%mult_type = 2
    !   mult(HL_S_LSL)%amat => mat(1:prim%groups_on_node,LSLrange)
    !   mult(HL_S_LSL)%bmat => mat(1:prim%groups_on_node,Srange)
    !   mult(HL_S_LSL)%cmat => mat(1:prim%groups_on_node,LHrange)
    !   mult(HL_S_LSL)%ahalo => halo(LSLrange)
    !   mult(HL_S_LSL)%chalo => halo(LHrange)
    !   mult(HL_S_LSL)%ltrans => ltrans(LSLrange)
    !   mult(HL_S_LSL)%bindex => Smatind
    !   mult(HL_S_LSL)%parts => parts
    !   mult(HL_S_LSL)%prim => prim
    !   mult(HL_S_LSL)%gcs => gcs
    !   call mult_ini(mult(HL_S_LSL),myid-1,prim%n_prim,parts)
    !endif
    ! 11-15
    mult(HLS_LS_L)%mult_type = 2
    mult(HLS_LS_L)%amat => mat(1:prim%groups_on_node,Lrange)
    mult(HLS_LS_L)%bmat => mat(1:prim%groups_on_node,SLrange)
    mult(HLS_LS_L)%cmat => mat(1:prim%groups_on_node,LSLrange)
    mult(HLS_LS_L)%ahalo => halo(Lrange)
    mult(HLS_LS_L)%chalo => halo(LSLrange)
    mult(HLS_LS_L)%ltrans => ltrans(Lrange)
    !mult(HLS_LS_L)%bindex => SLmatind
    mult(HLS_LS_L)%parts => parts
    mult(HLS_LS_L)%prim => prim
    mult(HLS_LS_L)%gcs => gcs
    call mult_ini(mult(HLS_LS_L),SLmatind,myid-1,prim%n_prim,parts)
    mult(SLS_LS_L)%mult_type = 2
    mult(SLS_LS_L)%amat => mat(1:prim%groups_on_node,Lrange)
    mult(SLS_LS_L)%bmat => mat(1:prim%groups_on_node,SLrange)
    mult(SLS_LS_L)%cmat => mat(1:prim%groups_on_node,SLSrange)
    mult(SLS_LS_L)%ahalo => halo(Lrange)
    mult(SLS_LS_L)%chalo => halo(SLSrange)
    mult(SLS_LS_L)%ltrans => ltrans(Lrange)
    !mult(SLS_LS_L)%bindex => SLmatind
    mult(SLS_LS_L)%parts => parts
    mult(SLS_LS_L)%prim => prim
    mult(SLS_LS_L)%gcs => gcs
    call mult_ini(mult(SLS_LS_L),SLmatind,myid-1,prim%n_prim,parts)
    mult(LHL_SL_S)%mult_type = 2
    mult(LHL_SL_S)%amat => mat(1:prim%groups_on_node,Srange)
    mult(LHL_SL_S)%bmat => mat(1:prim%groups_on_node,LSrange)
    mult(LHL_SL_S)%cmat => mat(1:prim%groups_on_node,SLSrange)
    mult(LHL_SL_S)%ahalo => halo(Srange)
    mult(LHL_SL_S)%chalo => halo(SLSrange)
    mult(LHL_SL_S)%ltrans => ltrans(Srange)
    !mult(LHL_SL_S)%bindex => LSmatind
    mult(LHL_SL_S)%parts => parts
    mult(LHL_SL_S)%prim => prim
    mult(LHL_SL_S)%gcs => gcs
    call mult_ini(mult(LHL_SL_S),LSmatind,myid-1,prim%n_prim,parts)
    mult(LSL_SL_H)%mult_type = 2
    mult(LSL_SL_H)%amat => mat(1:prim%groups_on_node,Hrange)
    mult(LSL_SL_H)%bmat => mat(1:prim%groups_on_node,LSrange)
    mult(LSL_SL_H)%cmat => mat(1:prim%groups_on_node,LSLrange)
    mult(LSL_SL_H)%ahalo => halo(Hrange)
    mult(LSL_SL_H)%chalo => halo(LSLrange)
    mult(LSL_SL_H)%ltrans => ltrans(Hrange)
    !mult(LSL_SL_H)%bindex => LSmatind
    mult(LSL_SL_H)%parts => parts
    mult(LSL_SL_H)%prim => prim
    mult(LSL_SL_H)%gcs => gcs
    call mult_ini(mult(LSL_SL_H),LSmatind,myid-1,prim%n_prim,parts)
    ra = rcut(SPrange)
    rc = rcut(Hrange)
    if(rc>=ra) then
       mult(SP_PS_H)%mult_type = 1
       mult(SP_PS_H)%amat => mat(1:prim%groups_on_node,SPrange)
       mult(SP_PS_H)%bmat => mat(1:prim%groups_on_node,PSrange)
       mult(SP_PS_H)%cmat => mat(1:prim%groups_on_node,Hrange)
       mult(SP_PS_H)%ahalo => halo(SPrange)
       mult(SP_PS_H)%chalo => halo(Hrange)
       mult(SP_PS_H)%ltrans => ltrans(SPrange)
       !mult(SP_PS_H)%bindex => PSmatind
       mult(SP_PS_H)%parts => parts
       mult(SP_PS_H)%prim => prim
       mult(SP_PS_H)%gcs => gcs
       call mult_ini(mult(SP_PS_H),PSmatind,myid-1,prim%n_prim,parts)
       mult(PAOP_PS_H)%mult_type = 1
       mult(PAOP_PS_H)%amat => mat(1:prim%groups_on_node,PAOPrange)
       mult(PAOP_PS_H)%bmat => mat(1:prim%groups_on_node,PSrange)
       mult(PAOP_PS_H)%cmat => mat(1:prim%groups_on_node,dHrange)
       mult(PAOP_PS_H)%ahalo => halo(PAOPrange)
       mult(PAOP_PS_H)%chalo => halo(dHrange)
       mult(PAOP_PS_H)%ltrans => ltrans(PAOPrange)
       !mult(PAOP_PS_H)%bindex => PSmatind
       mult(PAOP_PS_H)%parts => parts
       mult(PAOP_PS_H)%prim => prim
       mult(PAOP_PS_H)%gcs => gcs
       call mult_ini(mult(PAOP_PS_H),PSmatind,myid-1,prim%n_prim,parts)
    else
       mult(SP_PS_H)%mult_type = 2
       mult(SP_PS_H)%amat => mat(1:prim%groups_on_node,Hrange)
       mult(SP_PS_H)%bmat => mat(1:prim%groups_on_node,SPrange)
       mult(SP_PS_H)%cmat => mat(1:prim%groups_on_node,SPrange)
       mult(SP_PS_H)%ahalo => halo(Hrange)
       mult(SP_PS_H)%chalo => halo(SPrange)
       mult(SP_PS_H)%ltrans => ltrans(Hrange)
       !mult(SP_PS_H)%bindex => SPmatind
       mult(SP_PS_H)%parts => parts
       mult(SP_PS_H)%prim => prim
       mult(SP_PS_H)%gcs => gcs
       call mult_ini(mult(SP_PS_H),SPmatind,myid-1,prim%n_prim,parts)
       mult(PAOP_PS_H)%mult_type = 2
       mult(PAOP_PS_H)%amat => mat(1:prim%groups_on_node,dHrange)
       mult(PAOP_PS_H)%bmat => mat(1:prim%groups_on_node,SPrange)
       mult(PAOP_PS_H)%cmat => mat(1:prim%groups_on_node,PAOPrange)
       mult(PAOP_PS_H)%ahalo => halo(dHrange)
       mult(PAOP_PS_H)%chalo => halo(PAOPrange)
       mult(PAOP_PS_H)%ltrans => ltrans(dHrange)
       !mult(PAOP_PS_H)%bindex => SPmatind
       mult(PAOP_PS_H)%parts => parts
       mult(PAOP_PS_H)%prim => prim
       mult(PAOP_PS_H)%gcs => gcs
       call mult_ini(mult(PAOP_PS_H),SPmatind,myid-1,prim%n_prim,parts)
    endif
    ! 16-20
    mult(TS_T_T)%mult_type = 2
    mult(TS_T_T)%amat => mat(1:prim%groups_on_node,Trange)
    mult(TS_T_T)%bmat => mat(1:prim%groups_on_node,TTrrange)
    mult(TS_T_T)%cmat => mat(1:prim%groups_on_node,TSrange)
    mult(TS_T_T)%ahalo => halo(Trange)
    mult(TS_T_T)%chalo => halo(TSrange)
    mult(TS_T_T)%ltrans => ltrans(TTrrange)
    !mult(TS_T_T)%bindex => TTrmatind
    mult(TS_T_T)%parts => parts
    mult(TS_T_T)%prim => prim
    mult(TS_T_T)%gcs => gcs
    call mult_ini(mult(TS_T_T),TTrmatind,myid-1,prim%n_prim,parts)
    ra = rcut(THrange)
    rc = rcut(Lrange)
    if(rc>=ra) then
       mult(TH_T_L)%mult_type = 1
       mult(TH_T_L)%amat => mat(1:prim%groups_on_node,THrange)
       mult(TH_T_L)%bmat => mat(1:prim%groups_on_node,Trange)
       mult(TH_T_L)%cmat => mat(1:prim%groups_on_node,Lrange)
       mult(TH_T_L)%ahalo => halo(THrange)
       mult(TH_T_L)%chalo => halo(Lrange)
       mult(TH_T_L)%ltrans => ltrans(THrange)
       !mult(TH_T_L)%bindex => Tmatind
       mult(TH_T_L)%parts => parts
       mult(TH_T_L)%prim => prim
       mult(TH_T_L)%gcs => gcs
       call mult_ini(mult(TH_T_L),Tmatind,myid-1,prim%n_prim,parts)
    else
       mult(TH_T_L)%mult_type = 2
       mult(TH_T_L)%amat => mat(1:prim%groups_on_node,Lrange)
       mult(TH_T_L)%bmat => mat(1:prim%groups_on_node,TTrrange)
       mult(TH_T_L)%cmat => mat(1:prim%groups_on_node,THrange)
       mult(TH_T_L)%ahalo => halo(Lrange)
       mult(TH_T_L)%chalo => halo(THrange)
       mult(TH_T_L)%ltrans => ltrans(Lrange)
       !mult(TH_T_L)%bindex => TTrmatind
       mult(TH_T_L)%parts => parts
       mult(TH_T_L)%prim => prim
       mult(TH_T_L)%gcs => gcs
       call mult_ini(mult(TH_T_L),TTrmatind,myid-1,prim%n_prim,parts)
    endif
    mult(T_H_TH)%mult_type = 1
    mult(T_H_TH)%amat => mat(1:prim%groups_on_node,Trange)
    mult(T_H_TH)%bmat => mat(1:prim%groups_on_node,Hrange)
    mult(T_H_TH)%cmat => mat(1:prim%groups_on_node,THrange)
    mult(T_H_TH)%ahalo => halo(Trange)
    mult(T_H_TH)%chalo => halo(THrange)
    mult(T_H_TH)%ltrans => ltrans(Trange)
    !mult(T_H_TH)%bindex => Hmatind
    mult(T_H_TH)%parts => parts
    mult(T_H_TH)%prim => prim
    mult(T_H_TH)%gcs => gcs
    call mult_ini(mult(T_H_TH),Hmatind,myid-1,prim%n_prim,parts)
    mult(T_L_TL)%mult_type = 1
    mult(T_L_TL)%amat => mat(1:prim%groups_on_node,Trange)
    mult(T_L_TL)%bmat => mat(1:prim%groups_on_node,Lrange)
    mult(T_L_TL)%cmat => mat(1:prim%groups_on_node,TLrange)
    mult(T_L_TL)%ahalo => halo(Trange)
    mult(T_L_TL)%chalo => halo(TLrange)
    mult(T_L_TL)%ltrans => ltrans(Trange)
    !mult(T_L_TL)%bindex => Lmatind
    mult(T_L_TL)%parts => parts
    mult(T_L_TL)%prim => prim
    mult(T_L_TL)%gcs => gcs
    call mult_ini(mult(T_L_TL),Lmatind,myid-1,prim%n_prim,parts)
    mult(TL_T_L)%mult_type = 2
    mult(TL_T_L)%amat => mat(1:prim%groups_on_node,Lrange)
    mult(TL_T_L)%bmat => mat(1:prim%groups_on_node,TTrrange)
    mult(TL_T_L)%cmat => mat(1:prim%groups_on_node,TLrange)
    mult(TL_T_L)%ahalo => halo(Lrange)
    mult(TL_T_L)%chalo => halo(TLrange)
    mult(TL_T_L)%ltrans => ltrans(Lrange)
    !mult(TL_T_L)%bindex => TTrmatind
    mult(TL_T_L)%parts => parts
    mult(TL_T_L)%prim => prim
    mult(TL_T_L)%gcs => gcs
    call mult_ini(mult(TL_T_L),TTrmatind,myid-1,prim%n_prim,parts)
    !endif
!    call stop_timer(tmr_std_matrices)
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
!!  SOURCE
!!
  subroutine fmmi(prim)

    use basic_types
    use global_module, ONLY: iprint_mat
    use matrix_module, ONLY: deallocate_comms_data
    use matrix_data
    use maxima_module, ONLY: maxnabaprocs

    implicit none

    ! Passed variables
    type(primary_set) :: prim

    ! Local variables
    integer :: i,j

!    call start_timer(tmr_std_matrices)
   ! Multiplications
    nullify(mult(H_SP_SP)%amat,mult(H_SP_SP)%bmat,mult(H_SP_SP)%cmat,&
    mult(H_SP_SP)%ahalo,mult(H_SP_SP)%chalo,mult(H_SP_SP)%ltrans,&
    mult(H_SP_SP)%bindex,mult(H_SP_SP)%parts,mult(H_SP_SP)%prim,&
    mult(H_SP_SP)%gcs)
    nullify(mult(L_S_LS)%amat,mult(L_S_LS)%bmat,mult(L_S_LS)%cmat,&
    mult(L_S_LS)%ahalo,mult(L_S_LS)%chalo,mult(L_S_LS)%ltrans,&
    mult(L_S_LS)%bindex,mult(L_S_LS)%parts,mult(L_S_LS)%prim,&
    mult(L_S_LS)%gcs)
    nullify(mult(L_H_LH)%amat,mult(L_H_LH)%bmat,mult(L_H_LH)%cmat,&
    mult(L_H_LH)%ahalo,mult(L_H_LH)%chalo,mult(L_H_LH)%ltrans,&
    mult(L_H_LH)%bindex,mult(L_H_LH)%parts,mult(L_H_LH)%prim,&
    mult(L_H_LH)%gcs)
    nullify(mult(T_S_TS)%amat,mult(T_S_TS)%bmat,mult(T_S_TS)%cmat,&
    mult(T_S_TS)%ahalo,mult(T_S_TS)%chalo,mult(T_S_TS)%ltrans,&
    mult(T_S_TS)%bindex,mult(T_S_TS)%parts,mult(T_S_TS)%prim,&
    mult(T_S_TS)%gcs)
    !nullify(mult(S_TS_S)%amat,mult(S_TS_S)%bmat,mult(S_TS_S)%cmat,&
    !mult(S_TS_S)%ahalo,mult(S_TS_S)%chalo,mult(S_TS_S)%ltrans,&
    !mult(S_TS_S)%bindex,mult(S_TS_S)%parts,mult(S_TS_S)%prim,&
    !mult(S_TS_S)%gcs)
    !nullify(mult(S_LS_L)%amat,mult(S_LS_L)%bmat,mult(S_LS_L)%cmat,&
    !mult(S_LS_L)%ahalo,mult(S_LS_L)%chalo,mult(S_LS_L)%ltrans,&
    !mult(S_LS_L)%bindex,mult(S_LS_L)%parts,mult(S_LS_L)%prim,&
    !mult(S_LS_L)%gcs)
    nullify(mult(SL_S_SLS)%amat,mult(SL_S_SLS)%bmat,mult(SL_S_SLS)%cmat,&
    mult(SL_S_SLS)%ahalo,mult(SL_S_SLS)%chalo,mult(SL_S_SLS)%ltrans,&
    mult(SL_S_SLS)%bindex,mult(SL_S_SLS)%parts,mult(SL_S_SLS)%prim,&
    mult(SL_S_SLS)%gcs)
    nullify(mult(LH_L_SLS)%amat,mult(LH_L_SLS)%bmat,mult(LH_L_SLS)%cmat,&
    mult(LH_L_SLS)%ahalo,mult(LH_L_SLS)%chalo,mult(LH_L_SLS)%ltrans,&
    mult(LH_L_SLS)%bindex,mult(LH_L_SLS)%parts,mult(LH_L_SLS)%prim,&
    mult(LH_L_SLS)%gcs)
    nullify(mult(LS_L_LSL)%amat,mult(LS_L_LSL)%bmat,mult(LS_L_LSL)%cmat,&
    mult(LS_L_LSL)%ahalo,mult(LS_L_LSL)%chalo,mult(LS_L_LSL)%ltrans,&
    mult(LS_L_LSL)%bindex,mult(LS_L_LSL)%parts,mult(LS_L_LSL)%prim,&
    mult(LS_L_LSL)%gcs)
    nullify(mult(HL_S_LSL)%amat,mult(HL_S_LSL)%bmat,mult(HL_S_LSL)%cmat,&
    mult(HL_S_LSL)%ahalo,mult(HL_S_LSL)%chalo,mult(HL_S_LSL)%ltrans,&
    mult(HL_S_LSL)%bindex,mult(HL_S_LSL)%parts,mult(HL_S_LSL)%prim,&
    mult(HL_S_LSL)%gcs)
    nullify(mult(HLS_LS_L)%amat,mult(HLS_LS_L)%bmat,mult(HLS_LS_L)%cmat,&
    mult(HLS_LS_L)%ahalo,mult(HLS_LS_L)%chalo,mult(HLS_LS_L)%ltrans,&
    mult(HLS_LS_L)%bindex,mult(HLS_LS_L)%parts,mult(HLS_LS_L)%prim,&
    mult(HLS_LS_L)%gcs)
    !nullify(mult(LSL_SL_L)%amat,mult(LSL_SL_L)%bmat,mult(LSL_SL_L)%cmat,&
    !mult(LSL_SL_L)%ahalo,mult(LSL_SL_L)%chalo,mult(LSL_SL_L)%ltrans,&
    !mult(LSL_SL_L)%bindex,mult(LSL_SL_L)%parts,mult(LSL_SL_L)%prim,&
    !mult(LSL_SL_L)%gcs)
    nullify(mult(SLS_LS_L)%amat,mult(SLS_LS_L)%bmat,mult(SLS_LS_L)%cmat,&
    mult(SLS_LS_L)%ahalo,mult(SLS_LS_L)%chalo,mult(SLS_LS_L)%ltrans,&
    mult(SLS_LS_L)%bindex,mult(SLS_LS_L)%parts,mult(SLS_LS_L)%prim,&
    mult(SLS_LS_L)%gcs)
    nullify(mult(LHL_SL_S)%amat,mult(LHL_SL_S)%bmat,mult(LHL_SL_S)%cmat,&
    mult(LHL_SL_S)%ahalo,mult(LHL_SL_S)%chalo,mult(LHL_SL_S)%ltrans,&
    mult(LHL_SL_S)%bindex,mult(LHL_SL_S)%parts,mult(LHL_SL_S)%prim,&
    mult(LHL_SL_S)%gcs)
    nullify(mult(LSL_SL_H)%amat,mult(LSL_SL_H)%bmat,mult(LSL_SL_H)%cmat,&
    mult(LSL_SL_H)%ahalo,mult(LSL_SL_H)%chalo,mult(LSL_SL_H)%ltrans,&
    mult(LSL_SL_H)%bindex,mult(LSL_SL_H)%parts,mult(LSL_SL_H)%prim,&
    mult(LSL_SL_H)%gcs)
    nullify(mult(SP_PS_H)%amat,mult(SP_PS_H)%bmat,mult(SP_PS_H)%cmat,&
    mult(SP_PS_H)%ahalo,mult(SP_PS_H)%chalo,mult(SP_PS_H)%ltrans,&
    mult(SP_PS_H)%bindex,mult(SP_PS_H)%parts,mult(SP_PS_H)%prim,&
    mult(SP_PS_H)%gcs)
    nullify(mult(PAOP_PS_H)%amat,mult(PAOP_PS_H)%bmat,mult(PAOP_PS_H)%cmat,&
    mult(PAOP_PS_H)%ahalo,mult(PAOP_PS_H)%chalo,mult(PAOP_PS_H)%ltrans,&
    mult(PAOP_PS_H)%bindex,mult(PAOP_PS_H)%parts,mult(PAOP_PS_H)%prim,&
    mult(PAOP_PS_H)%gcs)
    nullify(mult(TS_T_T)%amat,mult(TS_T_T)%bmat,mult(TS_T_T)%cmat,&
    mult(TS_T_T)%ahalo,mult(TS_T_T)%chalo,mult(TS_T_T)%ltrans,&
    mult(TS_T_T)%bindex,mult(TS_T_T)%parts,mult(TS_T_T)%prim,&
    mult(TS_T_T)%gcs)
    nullify(mult(TH_T_L)%amat,mult(TH_T_L)%bmat,mult(TH_T_L)%cmat,&
    mult(TH_T_L)%ahalo,mult(TH_T_L)%chalo,mult(TH_T_L)%ltrans,&
    mult(TH_T_L)%bindex,mult(TH_T_L)%parts,mult(TH_T_L)%prim,&
    mult(TH_T_L)%gcs)
    nullify(mult(T_H_TH)%amat,mult(T_H_TH)%bmat,mult(T_H_TH)%cmat,&
    mult(T_H_TH)%ahalo,mult(T_H_TH)%chalo,mult(T_H_TH)%ltrans,&
    mult(T_H_TH)%bindex,mult(T_H_TH)%parts,mult(T_H_TH)%prim,&
    mult(T_H_TH)%gcs)
    nullify(mult(T_L_TL)%amat,mult(T_L_TL)%bmat,mult(T_L_TL)%cmat,&
    mult(T_L_TL)%ahalo,mult(T_L_TL)%chalo,mult(T_L_TL)%ltrans,&
    mult(T_L_TL)%bindex,mult(T_L_TL)%parts,mult(T_L_TL)%prim,&
    mult(T_L_TL)%gcs)
    nullify(mult(TL_T_L)%amat,mult(TL_T_L)%bmat,mult(TL_T_L)%cmat,&
    mult(TL_T_L)%ahalo,mult(TL_T_L)%chalo,mult(TL_T_L)%ltrans,&
    mult(TL_T_L)%bindex,mult(TL_T_L)%parts,mult(TL_T_L)%prim,&
    mult(TL_T_L)%gcs)
    ! Comms
    call deallocate_comms_data(mult(H_SP_SP)%comms)
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
    call deallocate_comms_data(mult(SP_PS_H)%comms)
    call deallocate_comms_data(mult(PAOP_PS_H)%comms)
    call deallocate_comms_data(mult(TS_T_T)%comms)
    call deallocate_comms_data(mult(TH_T_L)%comms)
    call deallocate_comms_data(mult(T_H_TH)%comms)
    call deallocate_comms_data(mult(T_L_TL)%comms)
    call deallocate_comms_data(mult(TL_T_L)%comms)
    do i=1,mx_trans
       do j = 1,gtrans(i)%n_rem_node!maxnabaprocs+1
          !if(associated(pairs(j,i)%submat)) deallocate(pairs(j,i)%submat)
          deallocate(pairs(j,i)%submat)
       end do
    end do
    deallocate(pairs)
    deallocate(Spairind,Lpairind,SPpairind,LSpairind,LHpairind,LSLpairind,Tpairind)
    !call dissociate_matrices
    ! Matrices
    call end_ops(prim,Srange,Smatind,S_trans)
    call end_ops(prim,dSrange,dSmatind)
    call end_ops(prim,Lrange,Lmatind,L_trans)
    call end_ops(prim,LTrrange,LTrmatind)
    call end_ops(prim,Hrange,Hmatind)
    call end_ops(prim,dHrange,dHmatind)
    call end_ops(prim,SPrange, SPmatind,SP_trans)
    call end_ops(prim,PAOPrange,PAOPmatind)
    call end_ops(prim,PSrange,PSmatind)
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
!    call stop_timer(tmr_std_matrices)
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
!!    Added timing calls and a LOCAL iprint_mult to activate them: we can make this global later
!!   2006/10/19 16:04 dave
!!    Changed iprint_mult to iprint_mat (global variable)
!!  SOURCE
!!
  subroutine LNV_matrix_multiply(electrons, energy, doK, doM1, doM2, doM3, doM4, dophi, doE, mat_M12, mat_M3, mat_M4, mat_phi)

    use numbers
    use datatypes
    use global_module, ONLY: iprint_mat, IPRINT_TIME_THRES3
    use matrix_data, ONLY: LHrange, SLSrange,LSLrange, Lrange, Srange, Hrange
    use GenComms, ONLY: gsum, inode, ionode, my_barrier, gmax,mtime
    use timer_module

    implicit none

    ! Passed variables
    real(double) :: electrons,energy
    integer :: myid
    integer :: mat_M12, mat_M3, mat_M4, mat_phi
    logical :: doK, doM1, doM2, doM3, doM4, dophi, doE

    ! This routine is basically calls to matrix operations,
    !   so these are timed within the routines themselves

    ! Local variables
    real(double) :: electrons_1, electrons_2, t0, t1
    integer :: matLH, matHL, matLHL, matHLS, matSLH, matA, matSLS, matLSL, matLHLSL, matLSLHL, matB, matC, matSLSLS, matLSLSL
    type(cq_timer) :: tmr_l_tmp1

    call start_timer(tmr_l_tmp1,WITH_LEVEL)
    if((doM1.OR.doM2).AND.mat_M12==0) call cq_abort("LNVmm needs proper matrix for M12")
    if(doM3.AND.mat_M3==0) call cq_abort("LNVmm needs proper matrix for M3")
    if(doM4.AND.mat_M4==0) call cq_abort("LNVmm needs proper matrix for M4")
    if(dophi.AND.mat_phi==0) call cq_abort("LNVmm needs proper matrix for phi")
     myid = inode-1
     if(iprint_mat>3) t0 = mtime()
     call matrix_product(matL, matS, matLS, mult(L_S_LS))
     if(iprint_mat>3) then
        t1 = mtime()
        if(inode==ionode) write(io_lun,*) 'LS time: ',t1-t0
        t0 = t1
     end if
     call matrix_transpose(matLS, matSL)
     if(iprint_mat>3) then
        t1 = mtime()
        if(inode==ionode) write(io_lun,*) 'LS trans time: ',t1-t0
        t0 = t1
     end if
     if(doM1.OR.doM2.OR.doM3) then
        matLH = allocate_temp_matrix(LHrange,LH_trans)
        matHL = allocate_temp_matrix(LHrange,LH_trans)
        if(iprint_mat>3) then
           t1 = mtime()
           if(inode==ionode) write(io_lun,*) 'LH alloc time: ',t1-t0
           t0 = t1
        end if
        call matrix_product(matL,matH,matLH,mult(L_H_LH))
        if(iprint_mat>3) then
           t1 = mtime()
           if(inode==ionode) write(io_lun,*) 'LH time: ',t1-t0
           t0 = t1
        end if
        call matrix_transpose(matLH, matHL)
        if(iprint_mat>3) then
           t1 = mtime()
           if(inode==ionode) write(io_lun,*) 'LH trans time: ',t1-t0
           t0 = t1
        end if
     endif
     if(doM1.OR.doM2) then ! This is type 1 or 2 depending on H
        matLHL = allocate_temp_matrix(SLSrange, 0)
        call matrix_product(matLH, matL, matLHL,mult(LH_L_SLS))
     endif
     if(doM2) then
        matLHLSL = allocate_temp_matrix(Srange, S_trans)
        matLSLHL = allocate_temp_matrix(Srange, S_trans)
        call matrix_product(matLHL, matLS, matLHLSL, mult(LHL_SL_S))
        call matrix_transpose(matLHLSL, matLSLHL)
     endif
     if(doM3) then
        if(iprint_mat>3) t0 = mtime()
        matHLS = allocate_temp_matrix(LSLrange,LSL_trans)
        matSLH = allocate_temp_matrix(LSLrange,LSL_trans)
        matA = allocate_temp_matrix(LSLrange,LSL_trans)
        matB = allocate_temp_matrix(Lrange,L_trans)
        matC = allocate_temp_matrix(Lrange,L_trans)
        if(iprint_mat>3) then
           t1 = mtime()
           if(inode==ionode) write(io_lun,*) 'Alloc time: ',t1-t0
           t0 = t1
        end if
        call matrix_product(matHL, matS, matHLS, mult(HL_S_LSL))
        if(iprint_mat>3) then
           t1 = mtime()
           if(inode==ionode) write(io_lun,*) 'HLS time: ',t1-t0
           t0 = t1
        end if
        call matrix_transpose(matHLS, matSLH)
        if(iprint_mat>3) then
           t1 = mtime()
           if(inode==ionode) write(io_lun,*) 'HLS trans time: ',t1-t0
           t0 = t1
        end if
        call matrix_sum(zero,matA,two,matHLS)
        call matrix_sum(one,matA,one,matSLH)
        if(iprint_mat>3) then
           t1 = mtime()
           if(inode==ionode) write(io_lun,*) 'A time: ',t1-t0
           t0 = t1
        end if
        call matrix_product(matA, matSL, matB, mult(HLS_LS_L))
        if(iprint_mat>3) then
           t1 = mtime()
           if(inode==ionode) write(io_lun,*) 'B time: ',t1-t0
           t0 = t1
        end if
        call matrix_transpose(matB, matC)
        if(iprint_mat>3) then
           t1 = mtime()
           if(inode==ionode) write(io_lun,*) 'B trans time: ',t1-t0
           t0 = t1
        end if
     endif
     if(dophi) then
        if(iprint_mat>3) t0 = mtime()
        matSLS = allocate_temp_matrix(SLSrange,0)
        matSLSLS = allocate_temp_matrix(Lrange,0)
        if(iprint_mat>3) then
           t1 = mtime()
           if(inode==ionode) write(io_lun,*) 'SLS alloc time: ',t1-t0
           t0 = t1
        end if
        call matrix_product(matSL,matS,matSLS,mult(SL_S_SLS))
        if(iprint_mat>3) then
           t1 = mtime()
           if(inode==ionode) write(io_lun,*) 'SLS time: ',t1-t0
           t0 = t1
        end if
        call matrix_product(matSLS,matSL,matSLSLS,mult(SLS_LS_L))
        if(iprint_mat>3) then
           t1 = mtime()
           if(inode==ionode) write(io_lun,*) 'SLSLS time: ',t1-t0
           t0 = t1
        end if
     endif
     if(doK.OR.doE.OR.doM4) then
        if(iprint_mat>3) t0 = mtime()
        matLSL = allocate_temp_matrix(LSLrange,0)
        matLSLSL = allocate_temp_matrix(Hrange,0)
        if(iprint_mat>3) then
           t1 = mtime()
           if(inode==ionode) write(io_lun,*) 'LSL alloc time: ',t1-t0
           t0 = t1
        end if
        call matrix_product(matLS,matL,matLSL,mult(LS_L_LSL))
        if(iprint_mat>3) then
           t1 = mtime()
           if(inode==ionode) write(io_lun,*) 'LSL time: ',t1-t0
           t0 = t1
        end if
        call matrix_product(matLSL,matLS,matLSLSL,mult(LSL_SL_H))
        if(iprint_mat>3) then
           t1 = mtime()
           if(inode==ionode) write(io_lun,*) 'LSLSL time: ',t1-t0
           t0 = t1
        end if
     endif
     if(doK.OR.doE) then
        call matrix_sum(zero,matK,three,matLSL)
        call matrix_sum(one,matK,-two,matLSLSL)
        if(iprint_mat>3) then
           t1 = mtime()
           if(inode==ionode) write(io_lun,*) 'K time: ',t1-t0
           t0 = t1
        end if
     endif
     if(doM4) then
        call matrix_sum(zero,mat_M4,24.0_double,matLSL)
        call matrix_sum(one,mat_M4,-24.0_double,matLSLSL)
     end if
     if(doE) then
        if(iprint_mat>3) t0 = mtime()
        energy = two*matrix_product_trace(matH,matK)
        if(iprint_mat>3) then
           t1 = mtime()
           if(inode==ionode) write(io_lun,*) 'energy dot time: ',t1-t0
           t0 = t1
        end if
     endif
     if(doM1.AND.doM2) then
        call matrix_sum(zero,mat_M12,three,matLHL)
        call matrix_sum(one,mat_M12,-two,matLHLSL)
        call matrix_sum(one,mat_M12,-two,matLSLHL)
     endif
     if(doM3) then
        if(iprint_mat>3) t0 = mtime()
        call matrix_sum(zero,mat_M3,six,matHLS)
        if(iprint_mat>3) then
           t1 = mtime()
           if(inode==ionode) write(io_lun,*) 'M3 time: ',t1-t0
           t0 = t1
        end if
        call matrix_sum(one,mat_M3,six,matSLH)
        if(iprint_mat>3) then
           t1 = mtime()
           if(inode==ionode) write(io_lun,*) 'M3 add time: ',t1-t0
           t0 = t1
        end if
        call matrix_sum(one,mat_M3,-two,matB)
        if(iprint_mat>3) then
           t1 = mtime()
           if(inode==ionode) write(io_lun,*) 'M3 axpy B time: ',t1-t0
           t0 = t1
        end if
        call matrix_sum(one,mat_M3,-two,matC)
        if(iprint_mat>3) then
           t1 = mtime()
           if(inode==ionode) write(io_lun,*) 'M3 axpy C time: ',t1-t0
           t0 = t1
        end if
     endif
     if(dophi) then
        if(iprint_mat>3) t0 = mtime()
        call matrix_sum(zero,mat_phi,one,matSLS)
        if(iprint_mat>3) then
           t1 = mtime()
           if(inode==ionode) write(io_lun,*) 'SLS prune time: ',t1-t0
           t0 = t1
        end if
        electrons_1 = matrix_product_trace(matL,mat_phi)
        electrons_2 = matrix_product_trace(matL,matSLSLS)
        electrons = 6.0_double * electrons_1 - four * electrons_2
        if(iprint_mat>3) then
           t1 = mtime()
           if(inode==ionode) write(io_lun,*) 'electrons dot time: ',t1-t0
           t0 = t1
        end if
        if(iprint_mat>3) then
           t1 = mtime()
           if(inode==ionode) write(io_lun,*) 'electrons gsum time: ',t1-t0
           t0 = t1
        end if
        call matrix_sum(12.0_double,mat_phi,-12.0_double,matSLSLS)
        if(iprint_mat>3) then
           t1 = mtime()
           if(inode==ionode) write(io_lun,*) 'phi time: ',t1-t0
           t0 = t1
        end if
     endif
     call my_barrier
     if(iprint_mat>3) t0 = mtime()
     if(doK.OR.doE.OR.doM4) then
        call free_temp_matrix(matLSLSL)
        call free_temp_matrix(matLSL)
     end if
     if(dophi) then
        call free_temp_matrix(matSLSLS)
        call free_temp_matrix(matSLS)
     end if
     if(doM3) then
        call free_temp_matrix(matC)
        call free_temp_matrix(matB)
        call free_temp_matrix(matA)
        call free_temp_matrix(matSLH)
        call free_temp_matrix(matHLS)
     end if
     if(doM2) then 
        call free_temp_matrix(matLSLHL)
        call free_temp_matrix(matLHLSL)
     end if
     if(doM1.OR.doM2) then 
        call free_temp_matrix(matLHL)
     end if
     if(doM1.OR.doM2.OR.doM3) then
        call free_temp_matrix(matHL)
        call free_temp_matrix(matLH)
     end if
     if(iprint_mat>3) then
        t1 = mtime()
        if(inode==ionode) write(io_lun,*) 'dealloc time: ',t1-t0
     end if
     call stop_print_timer(tmr_l_tmp1,"LNV_matrix_multiply",IPRINT_TIME_THRES3)
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
!!  SOURCE
!!
  subroutine symmetrise_L()

    use datatypes
    use numbers
    use matrix_data, ONLY: Lrange

    implicit none

    ! Local variables
    integer :: length
    integer :: mattmp

    mattmp = allocate_temp_matrix(Lrange, L_trans)
    call matrix_transpose(matL,mattmp)
    call matrix_sum(half,matL,half,mattmp)
    call free_temp_matrix(mattmp)
  end subroutine symmetrise_L
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
  subroutine end_ops(prim,range,matind, gtran)

    use basic_types
    use matrix_data, ONLY: mat, halo
    
    implicit none
    
    ! Passed variables
    type(primary_set) :: prim
    integer :: range
    integer, pointer, dimension(:) :: matind
    integer, OPTIONAL :: gtran
    ! Local variables
    integer :: i

    ! Global transpose
    if(PRESENT(gtran)) call deallocate_trans_rem(gtrans(gtran))
    ! Matrices
    call deallocate_trans(ltrans(range))
    call deallocate_halo(halo(range))
    call deallocate_matrix(mat(:,range),prim%groups_on_node)
    deallocate(matind)
    do i=1,prim%groups_on_node
       nullify(mat(i,range)%n_nab,mat(i,range)%i_seq,mat(i,range)%i_acc,&
            mat(i,range)%i_nd_acc,mat(i,range)%npxyz,mat(i,range)%ndimj)
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
!!   Associates the matrix pointers in private type with actual matrices in matrix_data
!!
!!   I realise that this is ugly; it does, however, mean that we can test the new matrix routines
!!   without too many changes
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
!!  SOURCE
!!
  subroutine associate_matrices

    use numbers, ONLY: zero
    use matrix_data, ONLY: Srange, Hrange, Lrange, Trange, PSrange,&
         LSrange, SPrange, SLrange, dSrange, dHrange, TTrrange, mat
    use global_module, ONLY: area_matrices, iprint_mat
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl
!%%!    use matrix_data, ONLY: data_S, data_H, data_L, data_T, data_Ttran, &
!%%!         data_LS, data_SL, data_SC, data_CS, data_K, data_M4, data_U, data_phi, &
!%%!         data_M12, data_UT, data_KE, data_NL, Srange, Hrange, Lrange, Trange, PSrange,&
!%%!         LSrange, SPrange, SLrange, dSrange, dHrange, data_dH, data_dS, mat,TTrrange
    use GenComms, ONLY: inode, ionode

    implicit none

    integer :: stat, i

    ! Assign indices for the matrices that we need
    matS     = 1
    matH     = 2
    matL     = 3
    matT     = 4
    matTtran = 5
    matLS    = 6
    matSL    = 7
    matSC    = 8
    matCS    = 9
    matK     = 10
    matM4    = 11
    matU     = 12
    matphi   = 13
    matM12   = 14
    matUT    = 15
    matKE    = 16
    matNL    = 17
    matdH    = 18
    matdS    = 19
    current_matrix = 19
!%%!    ! Now associate pointers with correct arrays
!%%!    mat_p(matS  )%matrix => data_S
!%%!    mat_p(matH  )%matrix => data_H
!%%!    mat_p(matL  )%matrix => data_L
!%%!    mat_p(matT  )%matrix => data_T
!%%!    mat_p(matTtran)%matrix => data_Ttran
!%%!    mat_p(matLS )%matrix => data_LS
!%%!    mat_p(matSL )%matrix => data_SL
!%%!    mat_p(matSC )%matrix => data_SC
!%%!    mat_p(matCS )%matrix => data_CS
!%%!    mat_p(matK  )%matrix => data_K
!%%!    mat_p(matM4 )%matrix => data_M4
!%%!    mat_p(matU  )%matrix => data_U
!%%!    mat_p(matphi)%matrix => data_phi
!%%!    mat_p(matM12)%matrix => data_M12
!%%!    mat_p(matUT )%matrix => data_UT
!%%!    mat_p(matKE )%matrix => data_KE
!%%!    mat_p(matNL )%matrix => data_NL
!%%!    mat_p(matdH  )%matrix => data_dH
!%%!    mat_p(matdS  )%matrix => data_dS
    ! Set the dimensions of the arrays
    ! Type for f1
    mat_p(matS  )%sf1_type = sf
    mat_p(matH  )%sf1_type = sf
    mat_p(matL  )%sf1_type = sf
    mat_p(matT  )%sf1_type = sf
    mat_p(matTtran)%sf1_type = sf
    mat_p(matLS )%sf1_type = sf
    mat_p(matSL )%sf1_type = sf
    mat_p(matSC )%sf1_type = sf
    mat_p(matCS )%sf1_type = nlpf
    mat_p(matK  )%sf1_type = sf
    mat_p(matM4 )%sf1_type = sf
    mat_p(matU  )%sf1_type = sf
    mat_p(matphi)%sf1_type = sf
    mat_p(matM12)%sf1_type = sf
    mat_p(matUT )%sf1_type = nlpf
    mat_p(matKE )%sf1_type = sf
    mat_p(matNL )%sf1_type = sf
    mat_p(matdS )%sf1_type = paof
    mat_p(matdH )%sf1_type = paof
    ! Type for f2
    mat_p(matS  )%sf2_type = sf
    mat_p(matH  )%sf2_type = sf
    mat_p(matL  )%sf2_type = sf
    mat_p(matT  )%sf2_type = sf
    mat_p(matTtran)%sf2_type = sf
    mat_p(matLS )%sf2_type = sf
    mat_p(matSL )%sf2_type = sf
    mat_p(matSC )%sf2_type = nlpf
    mat_p(matCS )%sf2_type = sf
    mat_p(matK  )%sf2_type = sf
    mat_p(matM4 )%sf2_type = sf
    mat_p(matU  )%sf2_type = nlpf
    mat_p(matphi)%sf2_type = sf
    mat_p(matM12)%sf2_type = sf
    mat_p(matUT )%sf2_type = sf
    mat_p(matKE )%sf2_type = sf
    mat_p(matNL )%sf2_type = sf
    mat_p(matdS )%sf2_type = sf
    mat_p(matdH )%sf2_type = sf
    ! Set up index translation for ranges
    matrix_index(matS  ) = Srange  
    matrix_index(matH  ) = Hrange  
    matrix_index(matL  ) = Lrange  
    matrix_index(matT  ) = Trange  
    matrix_index(matTtran) = TTrrange  
    matrix_index(matLS ) = LSrange 
    matrix_index(matSL ) = SLrange 
    matrix_index(matSC ) = SPrange 
    matrix_index(matCS ) = PSrange 
    matrix_index(matK  ) = Hrange  
    matrix_index(matM4 ) = Srange  
    matrix_index(matU  ) = SPrange 
    matrix_index(matphi) = Lrange  
    matrix_index(matM12) = Srange  
    matrix_index(matUT ) = PSrange 
    matrix_index(matKE ) = Hrange  
    matrix_index(matNL ) = Hrange  
    matrix_index(matdS  ) = dSrange 
    matrix_index(matdH  ) = dHrange 
    ! Real lengths
    ! Length
!%%!    mat_p(matS  )%length = nsf*nsf*mx_at_prim*mx_snab
!%%!    mat_p(matH  )%length = nsf*nsf*mx_at_prim*mx_hnab
!%%!    mat_p(matL  )%length = nsf*nsf*mx_at_prim*mx_lnab
!%%!    mat_p(matT  )%length = nsf*nsf*mx_at_prim*mx_tnab
!%%!    mat_p(matTtran)%length = nsf*nsf*mx_at_prim*mx_tnab
!%%!    mat_p(matLS )%length = nsf*nsf*mx_at_prim*mx_lsnab
!%%!    mat_p(matSL )%length = nsf*nsf*mx_at_prim*mx_lsnab
!%%!    mat_p(matSC )%length = nsf*ncf*mx_at_prim*mx_scnab
!%%!    mat_p(matCS )%length = ncf*nsf*mx_at_prim*mx_scnab
!%%!    mat_p(matK  )%length = nsf*nsf*mx_at_prim*mx_hnab
!%%!    mat_p(matM4 )%length = nsf*nsf*mx_at_prim*mx_snab
!%%!    mat_p(matU  )%length = nsf*ncf*mx_at_prim*mx_scnab
!%%!    mat_p(matphi)%length = nsf*nsf*mx_at_prim*mx_lnab
!%%!    mat_p(matM12)%length = nsf*nsf*mx_at_prim*mx_snab
!%%!    mat_p(matUT )%length = ncf*nsf*mx_at_prim*mx_scnab
!%%!    mat_p(matKE )%length = nsf*nsf*mx_at_prim*mx_hnab
!%%!    mat_p(matNL )%length = nsf*nsf*mx_at_prim*mx_hnab
!%%!    mat_p(matdS  )%length = npao*nsf*mx_at_prim*mx_snab
!%%!    mat_p(matdH  )%length = npao*nsf*mx_at_prim*mx_hnab
    do i=1,current_matrix
       if(mat_p(i)%length==0) then
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
          if(stat/=0) call cq_abort("Error allocating matrix ",i,mat_p(i)%length)
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
    if(iprint_mat>3.AND.inode==ionode) then
       do stat=1,current_matrix
          write(io_lun,fmt='(2x,"Proc: ",i5," Matrix no., length: ",i4,i8)') stat,mat_p(stat)%length
       end do
    end if
    ! Set up index translation for transposes
    trans_index(matS  ) = S_trans 
    trans_index(matH  ) = 0
    trans_index(matL  ) = L_trans 
    trans_index(matT  ) = T_trans 
    trans_index(matTtran) = T_trans 
    trans_index(matLS ) = LS_trans
    trans_index(matSL ) = LS_trans
    trans_index(matSC ) = SP_trans
    trans_index(matCS ) = SP_trans
    trans_index(matK  ) = 0
    trans_index(matM4 ) = S_trans 
    trans_index(matU  ) = SP_trans
    trans_index(matphi) = L_trans 
    trans_index(matM12) = S_trans 
    trans_index(matUT ) = SP_trans
    trans_index(matKE ) = 0
    trans_index(matNL ) = 0
    trans_index(matdS  ) = 0
    trans_index(matdH  ) = 0
  end subroutine associate_matrices
!!***

  subroutine dissociate_matrices

    implicit none

    integer :: stat

    call start_timer(tmr_std_allocation)
    deallocate(mat_p(matdS )%matrix,STAT=stat) !19
    if(stat/=0) call cq_abort("Error deallocating matrix dH ",mat_p(matdS)%length)
    deallocate(mat_p(matdH )%matrix,STAT=stat) !18
    if(stat/=0) call cq_abort("Error deallocating matrix dS ",mat_p(matdH)%length)
    deallocate(mat_p(matNL )%matrix,STAT=stat) !17
    if(stat/=0) call cq_abort("Error deallocating matrix NL ",mat_p(matNL)%length)
    deallocate(mat_p(matKE )%matrix,STAT=stat) !16
    if(stat/=0) call cq_abort("Error deallocating matrix KE ",mat_p(matKE)%length)
    deallocate(mat_p(matUT )%matrix,STAT=stat) !15
    if(stat/=0) call cq_abort("Error deallocating matrix UT ",mat_p(matUT)%length)
    deallocate(mat_p(matM12)%matrix,STAT=stat) !14
    if(stat/=0) call cq_abort("Error deallocating matrix M12 ",mat_p(matM12)%length)
    deallocate(mat_p(matphi)%matrix,STAT=stat) !13
    if(stat/=0) call cq_abort("Error deallocating matrix phi ",mat_p(matphi)%length)
    deallocate(mat_p(matU  )%matrix,STAT=stat) !12
    if(stat/=0) call cq_abort("Error deallocating matrix U ",mat_p(matU)%length)
    deallocate(mat_p(matM4 )%matrix,STAT=stat) !11
    if(stat/=0) call cq_abort("Error deallocating matrix M4 ",mat_p(matM4)%length)
    deallocate(mat_p(matK  )%matrix,STAT=stat) !10
    if(stat/=0) call cq_abort("Error deallocating matrix K ",mat_p(matK)%length)
    deallocate(mat_p(matCS )%matrix,STAT=stat) !9
    if(stat/=0) call cq_abort("Error deallocating matrix CS ",mat_p(matCS)%length)
    deallocate(mat_p(matSC )%matrix,STAT=stat) !8
    if(stat/=0) call cq_abort("Error deallocating matrix SC ",mat_p(matSC)%length)
    deallocate(mat_p(matSL )%matrix,STAT=stat) !7
    if(stat/=0) call cq_abort("Error deallocating matrix SL ",mat_p(matSL)%length)
    deallocate(mat_p(matLS )%matrix,STAT=stat) !6
    if(stat/=0) call cq_abort("Error deallocating matrix LS ",mat_p(matLS)%length)
    deallocate(mat_p(matTtran)%matrix,STAT=stat) !5
    if(stat/=0) call cq_abort("Error deallocating matrix Ttran ",mat_p(matTtran)%length)
    deallocate(mat_p(matT  )%matrix,STAT=stat) !4
    if(stat/=0) call cq_abort("Error deallocating matrix T ",mat_p(matT)%length)
    deallocate(mat_p(matL  )%matrix,STAT=stat) !3
    if(stat/=0) call cq_abort("Error deallocating matrix L ",mat_p(matL)%length)
    deallocate(mat_p(matH  )%matrix,STAT=stat) !2
    if(stat/=0) call cq_abort("Error deallocating matrix H ",mat_p(matH)%length)
    deallocate(mat_p(matS  )%matrix,STAT=stat) !1
    if(stat/=0) call cq_abort("Error deallocating matrix S ",mat_p(matS)%length)
    call stop_timer(tmr_std_allocation)
  end subroutine dissociate_matrices

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
    use matrix_module, ONLY: matrix_mult
    use multiply_module, ONLY: mat_mult
    use GenComms, ONLY: myid

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
    if(PRESENT(debug)) then
       !write(io_lun,*) 'Debug set in matrix_product'
       if(mult%mult_type==1) then
          call mat_mult(myid,mat_p(A)%matrix,mat_p(A)%length,mat_p(B)%matrix,mat_p(B)%length, &
               mat_p(C)%matrix,mat_p(C)%length,mult,debug)
       else
          call mat_mult(myid,mat_p(C)%matrix,mat_p(C)%length,mat_p(B)%matrix,mat_p(B)%length, &
               mat_p(A)%matrix,mat_p(A)%length,mult,debug)
       endif
     else
        if(mult%mult_type==1) then
          call mat_mult(myid,mat_p(A)%matrix,mat_p(A)%length,mat_p(B)%matrix,mat_p(B)%length, &
               mat_p(C)%matrix,mat_p(C)%length,mult)
       else
          call mat_mult(myid,mat_p(C)%matrix,mat_p(C)%length,mat_p(B)%matrix,mat_p(B)%length, &
               mat_p(A)%matrix,mat_p(A)%length,mult)
       endif
    end if
!    call stop_timer(tmr_std_matrices)
  end subroutine matrix_product
!!***

  subroutine matrix_transpose(A, AT)

    use mpi
    use GenComms, ONLY: myid
    use matrix_data, ONLY: halo
    use global_module, ONLY: sf
    use maxima_module, ONLY: maxatomsproc, maxnabaprocs
    use multiply_module, ONLY: loc_trans, mat_tran
    use comms_module, ONLY: send_trans_data

    implicit none

    ! Passed variables
    integer :: A, AT
    ! Local variables
    integer :: mat_temp
    integer  :: nreqs(1:1000), isend, ii, ierr=0
    integer, dimension(MPI_STATUS_SIZE) :: mpi_stat

!    call start_timer(tmr_std_matrices)
    ! Allocate temporary matrix to receive local transpose
    if(mat_p(A)%sf2_type/=sf.OR.mat_p(A)%sf1_type/=sf) then
       mat_temp = allocate_temp_matrix(matrix_index(A),0,mat_p(A)%sf2_type,mat_p(A)%sf1_type)
    else
       mat_temp = allocate_temp_matrix(matrix_index(A),0)
    end if
    ! Do local transpose for A
    !write(io_lun,*) myid,' Doing loc_trans'
    call loc_trans( ltrans(matrix_index(A)), halo(matrix_index(A)),mat_p(A)%matrix,mat_p(A)%length, &
         mat_p(mat_temp)%matrix,mat_p(mat_temp)%length,0)
    ! Send out the data
    !write(io_lun,*) myid,' Doing send_trans'
    call send_trans_data(maxatomsproc,gtrans(trans_index(A)),mat_p(mat_temp)%matrix,&
         myid+1, isend, nreqs,pairs(:,trans_index(A)),maxnabaprocs+1)
    ! Receive it and built AT
    !write(io_lun,*) myid,' Doing mat_tran'
    call mat_tran(maxatomsproc,myid+1,gtrans(trans_index(A)), &
         mat_p(AT)%matrix,mat_p(mat_temp)%matrix,pairs(:,trans_index(A)),maxnabaprocs+1)

    !MPI_Wait  T. Miyazaki 2003/11/18
    do ii= 1, isend
       call MPI_Wait(nreqs(ii),mpi_stat,ierr)
       if(ierr/=0) call cq_abort("Error in mat_trans: waiting for send to finish",ii)
    enddo
    call free_temp_matrix(mat_temp)
!    call stop_timer(tmr_std_matrices)
    return

  end subroutine matrix_transpose

  !!   This is equivalent to axpy: A = alpha*A + beta*B
  subroutine matrix_sum(alpha,A,beta,B)

    use datatypes
    use multiply_module, ONLY: matrix_add
    use matrix_data, ONLY: mat
    use primary_module, ONLY: bundle
    use GenComms, ONLY: myid
    use GenBlas, ONLY: axpy, scal
    
    implicit none

    ! Passed variables
    integer :: A, B
    real(double) :: alpha, beta

!    call start_timer(tmr_std_matrices)
    !write(io_lun,*) 'Summing matrices, ranges: ',A,B,matrix_index(A),matrix_index(B)
    if(matrix_index(A)==matrix_index(B)) then ! The matrices are the same range
       !write(io_lun,*) myid,' Scaling ',size(mat_p(A)%matrix),mat_p(A)%length
       call scal(mat_p(A)%length,alpha,mat_p(A)%matrix,1)
       !mat_p(A)%matrix = alpha*mat_p(A)%matrix 
       !write(io_lun,*) myid,' Adding ',size(mat_p(A)%matrix),size(mat_p(B)%matrix),mat_p(A)%length,mat_p(B)%length
       call axpy(mat_p(A)%length,beta,mat_p(B)%matrix,1,mat_p(A)%matrix,1)
    else
       if(matrix_index(A)<1.OR.matrix_index(A)>mx_matrices) call cq_abort("Matrix error in matrix_sum: ",matrix_index(A),A)
       if(matrix_index(B)<1.OR.matrix_index(B)>mx_matrices) call cq_abort("Matrix error in matrix_sum: ",matrix_index(B),B)
       call scal(mat_p(A)%length,alpha,mat_p(A)%matrix,1)
       call matrix_add(alpha,mat_p(A)%matrix,mat_p(A)%length,mat(:,matrix_index(A)), &
            beta,mat_p(B)%matrix,mat_p(B)%length,mat(:,matrix_index(B)),bundle,myid)
    end if
!    call stop_timer(tmr_std_matrices)
    return
  end subroutine matrix_sum

  subroutine matrix_scale(alpha,A)

    use datatypes
    use GenBlas, ONLY: scal

    implicit none

    integer :: A
    real(double) :: alpha

!    call start_timer(tmr_std_matrices)
    ! We could do this with a call to scal (the BLAS routine)
    if(A<1.OR.A>current_matrix) call cq_abort("Error in matrix_scale: ",A,mx_matrices)
    !write(io_lun,*) 'Scaling: ',A,alpha,size(mat_p(A)%matrix)
    call scal(mat_p(A)%length,alpha,mat_p(A)%matrix,1)
    !mat_p(A)%matrix = alpha*mat_p(A)%matrix
!    call stop_timer(tmr_std_matrices)

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
!!   Note: in writing this, we have assumed that we'll only ever take the trace of a matrix
!!   which is <sf|A|sf> (i.e. that it has support functions on-site).  This is explicit in the
!!   use of nsf_species for the inner loop.  It would be perfectly possible to change this with
!!   a switch based on the type of matrix if that became necessary
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
!!  SOURCE
!!
  real(double) function matrix_trace(A)

    use datatypes
    use numbers
    use matrix_data, ONLY: mat
    use primary_module, ONLY: bundle
    use GenComms, ONLY: gsum

    implicit none

    ! Passed variables
    integer :: A

    ! Local variables
    integer :: np,i,j, ierr, iprim, sf1, loc, Ah

!    call start_timer(tmr_std_matrices)
    Ah = matrix_index(A)
    matrix_trace = zero
    iprim = 0
    do np = 1,bundle%groups_on_node  ! Loop over partitions on processor
       if(bundle%nm_nodgroup(np)>0) then
          do i=1,bundle%nm_nodgroup(np) ! Loop over atoms in partition
             iprim = iprim + 1
             sf1 = mat(np,Ah)%ndimi(i)
             loc = (sf-1)*sf1+sf-1+mat(np,Ah)%onsite(i)
             if(loc>mat_p(A)%length) call cq_abort("Overflow error in trace: ",loc,mat_p(A)%length)
             do j=1,sf1
                loc = (j-1)*sf1+j-1+mat(np,Ah)%onsite(i)
                matrix_trace = matrix_trace + mat_p(A)%matrix(loc) ! Add
             enddo  ! End loop over NSF
          enddo ! End loop over atoms in partition
       endif ! End if atoms in partition
    enddo ! End loop over partitions on processor
    call gsum(matrix_trace) ! Global sum
!    call stop_timer(tmr_std_matrices)
  end function matrix_trace
!!***

!! Atom-by-atom matrix_product_trace e.g. for atom charges
!!   Call as atom_trace(A,B,prim%n_prim,res) where A and B are matrix labels (e.g. matK, matS),
!!   len is bundle%n_prim and res is a real(double) array of size len passed in.
!! D.R. Bowler, 2008/03/21
!!
  subroutine atom_trace(A,B,len,res)

    use datatypes
    use numbers, ONLY: zero
    use matrix_data, ONLY: mat
    use primary_module, ONLY: bundle
    use GenBlas, ONLY: dot

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
    if(matrix_index(A)==matrix_index(B)) then           ! The matrices are the same range
       do np = 1,bundle%groups_on_node                  ! Loop over partitions on processor
          if(bundle%nm_nodgroup(np)>0) then              ! if there are atoms in partition
             do i=1,bundle%nm_nodgroup(np)                ! loop over atoms in partition
                iprim = iprim + 1   ! count no. of atoms
                ist = mat(np,Ah)%nd_offset + mat(np,Ah)%i_nd_acc(i)
                if(i<bundle%nm_nodgroup(np)) then         ! if still going through partition atoms:
                   l = mat(np,Ah)%i_nd_acc(i+1) - mat(np,Ah)%i_nd_acc(i)
                else                                      ! if finished all atoms in partition
                   l = mat(np,Ah)%part_nd_nabs - mat(np,Ah)%i_nd_acc(i) + 1
                end if
                ! atomic charge for iprim atom:
                res(iprim) = dot(l,mat_p(A)%matrix(ist:ist+l-1),1,mat_p(B)%matrix(ist:ist+l-1),1)
             end do
          end if
       end do ! recovered the charges per processor, need to gsum
    end if
!    call stop_timer(tmr_std_matrices)

  end subroutine atom_trace
!!***


  real(double) function matrix_product_trace(A,B)

    use datatypes
    use numbers
    use matrix_data, ONLY: mat
    use GenBlas, ONLY: dot, vdot, asum

    implicit none

    ! Passed variables
    integer :: A,B

    ! Local variables
    integer :: np,i,j, ierr

!    call start_timer(tmr_std_matrices)
    if(matrix_index(A)==matrix_index(B)) then ! The matrices are the same range
       if(mat_p(A)%length/=mat_p(B)%length) call cq_abort("Length error in mat_prod_tr: ",mat_p(A)%length,mat_p(B)%length)
       matrix_product_trace = vdot(mat_p(A)%length,mat_p(A)%matrix,1,mat_p(B)%matrix,1)
    else ! For now I'm going to abort; we could have (say) optional argument giving the mult-type and range of the product
       call cq_abort("Matrix product for two different ranges not implemented")
    endif
!    call stop_timer(tmr_std_matrices)
  end function matrix_product_trace

  integer function return_matrix_len(A)

    implicit none

    ! Passed variable
    integer :: A
    return_matrix_len = mat_p(A)%length
  end function return_matrix_len

  ! i is the primary set number
  ! j is gcs%icover_ibeg(np)+ni - 1
  real(double) function return_matrix_value(A,np,ni,ip,nabj,f1,f2, onsite)

    use datatypes
    use basic_types, ONLY: primary_set, cover_set
    use matrix_module, ONLY: matrix_halo, matrix
    use matrix_data, ONLY: halo, mat
    use cover_module, ONLY: BCS_parts

    implicit none

    ! Passed variables
    integer :: A, np, ni, ip, nabj, f1, f2 
    integer, OPTIONAL :: onsite
    ! Local variables
    integer :: Ah, pos, jnp, jni, j, sf1

!    call start_timer(tmr_std_matrices)
    Ah = matrix_index(A)
    sf1 = mat(np,Ah)%ndimi(ni)
    if(PRESENT(onsite)) then
       pos = mat(np,Ah)%onsite(ni)
    else
       jnp = mat(np,Ah)%i_part(mat(np,Ah)%i_acc(ni)+nabj-1)
       jni = mat(np,Ah)%i_seq(mat(np,Ah)%i_acc(ni)+nabj-1)
       j = BCS_parts%icover_ibeg(jnp)+jni-1
       pos = halo(Ah)%i_h2d( (ip-1)*halo(Ah)%ni_in_halo + halo(Ah)%i_halo(j) )
    end if
    ! In future we'll just have matrix_pos here
    if((f2-1)*sf1+f1-1+pos>mat_p(A)%length) call cq_abort("Error in return_matrix_value",(f2-1)*sf1+f1-1+pos,mat_p(A)%length)
    return_matrix_value = mat_p(A)%matrix((f2-1)*sf1+f1-1+pos)
!    call stop_timer(tmr_std_matrices)
    !if(pos/=posin) write(io_lun,*) 'Pos: ',pos,posin,i,j,halo(Ah)%i_halo(j),Ah
  end function return_matrix_value

  real(double) function return_matrix_value_pos(A,i)

    use datatypes
    use GenComms, ONLY: cq_abort

    implicit none

    ! Passed variables
    integer, intent(in) :: A, i

!    call start_timer(tmr_std_matrices)
    if(i>mat_p(A)%length) call cq_abort("Matrix overrun in return_matrix_value_pos: ",i,A)
    ! In future we'll just have matrix_pos here
    return_matrix_value_pos = mat_p(A)%matrix(i)
  end function return_matrix_value_pos

  subroutine return_matrix_block_pos(A,i,val,size)

    use datatypes
    use GenComms, ONLY: cq_abort

    implicit none

    ! Passed variables
    integer :: A, i,size
    real(double), dimension(size) :: val

!    call start_timer(tmr_std_matrices)
    if(i+size-1>mat_p(A)%length) call cq_abort("Matrix overrun in return_matrix_value_pos: ",i,A)
    ! In future we'll just have matrix_pos here
    val(1:size) = mat_p(A)%matrix(i:i+size-1)
!    call stop_timer(tmr_std_matrices)
  end subroutine return_matrix_block_pos

  integer function matrix_pos(A,iprim,j_in_halo,f1,f2)

    use matrix_module, ONLY: matrix_halo
    use matrix_data, ONLY: halo
    use GenComms, ONLY: cq_abort

    implicit none

    ! Passed variables
    integer, intent(in) :: A, iprim, j_in_halo
    integer, OPTIONAL, intent(in) :: f1, f2
    ! Local variables
    integer :: Ah, sf1

!    call start_timer(tmr_std_matrices)
    if(A<1) call cq_abort("Error in matrix_pos: matrix number illegal: ",A)
    Ah = matrix_index(A)

    if(j_in_halo>halo(Ah)%mx_halo) call cq_abort('Halo indexing error in matrix_pos: ',j_in_halo,halo(Ah)%mx_halo)
    if(present(f1).AND.present(f2)) then
       sf1 = halo(Ah)%ndimi(iprim)
       matrix_pos = (f2-1)*sf1+f1-1+halo(Ah)%i_h2d( (iprim-1)*halo(Ah)%ni_in_halo + j_in_halo )
    else
       matrix_pos = halo(Ah)%i_h2d( (iprim-1)*halo(Ah)%ni_in_halo + j_in_halo )
    end if
    ! Some length check here ?
    if(matrix_pos>mat_p(A)%length) call cq_abort("Overrun error in matrix_pos: ",matrix_pos,mat_p(A)%length)
!    call stop_timer(tmr_std_matrices)
  end function matrix_pos

  subroutine scale_matrix_value(A,np,ni,ip,nabj,f1,f2,val,onsite)

    use datatypes
    use matrix_module, ONLY: matrix_halo, matrix
    use matrix_data, ONLY: halo, mat
    use cover_module, ONLY: BCS_parts

    implicit none

    ! Passed variables
    integer :: A, np, ni, ip, nabj, f1, f2
    integer, OPTIONAL :: onsite
    real(double) :: val
    ! Local variables
    integer :: Ah, pos, j, jnp, jni, sf1

!    call start_timer(tmr_std_matrices)
    Ah = matrix_index(A)
    sf1 = halo(Ah)%ndimi(ip)
    if(PRESENT(onsite)) then
       pos = mat(np,Ah)%onsite(ni)
    else
       jnp = mat(np,Ah)%i_part(mat(np,Ah)%i_acc(ni)+nabj-1)
       jni = mat(np,Ah)%i_seq(mat(np,Ah)%i_acc(ni)+nabj-1)
       j = BCS_parts%icover_ibeg(jnp)+jni-1
       pos = halo(Ah)%i_h2d( (ip-1)*halo(Ah)%ni_in_halo + halo(Ah)%i_halo(j) )
    end if
    ! In future we'll just have matrix_pos here
    if(f1-1+(f2-1)*sf1+pos>mat_p(A)%length) call cq_abort("Overrun error in scale_matrix_value: ",&
         f1-1+(f2-1)*sf1+pos,mat_p(A)%length)
    mat_p(A)%matrix(f1-1+(f2-1)*sf1+pos) = mat_p(A)%matrix(f1-1+(f2-1)*sf1+pos)*val
!    call stop_timer(tmr_std_matrices)
  end subroutine scale_matrix_value

  subroutine store_matrix_value(A,np,ni,ip,nabj,f1,f2,val,onsite)

    use datatypes
    use matrix_module, ONLY: matrix_halo, matrix
    use matrix_data, ONLY: halo, mat
    use cover_module, ONLY: BCS_parts
    use GenComms, ONLY: cq_abort

    implicit none

    ! Passed variables
    integer :: A, np, ni, ip, nabj, f1, f2
    integer, OPTIONAL :: onsite
    real(double) :: val
    ! Local variables
    integer :: Ah, pos, j, jnp, jni, sf1

!    call start_timer(tmr_std_matrices)
    Ah = matrix_index(A)
    sf1 = halo(Ah)%ndimi(ip)
    if(PRESENT(onsite)) then
       pos = mat(np,Ah)%onsite(ni)
    else
       jnp = mat(np,Ah)%i_part(mat(np,Ah)%i_acc(ni)+nabj-1)
       jni = mat(np,Ah)%i_seq(mat(np,Ah)%i_acc(ni)+nabj-1)
       j = BCS_parts%icover_ibeg(jnp)+jni-1
       pos = halo(Ah)%i_h2d( (ip-1)*halo(Ah)%ni_in_halo + halo(Ah)%i_halo(j) )
    end if
    if(f1-1+(f2-1)*sf1+pos>mat_p(A)%length) call cq_abort("Overrun error in store_matrix_value: ",&
         f1-1+(f2-1)*sf1+pos,mat_p(A)%length)
    ! In future we'll just have matrix_pos here
    mat_p(A)%matrix(f1-1+(f2-1)*sf1+pos) = mat_p(A)%matrix(f1-1+(f2-1)*sf1+pos) + val
!    call stop_timer(tmr_std_matrices)
  end subroutine store_matrix_value

  subroutine store_matrix_value_pos(A,i,val)

    use datatypes
    use GenComms, ONLY: cq_abort

    implicit none

    ! Passed variables
    integer, intent(in) :: A, i
    real(double), intent(in) :: val

!    call start_timer(tmr_std_matrices)
    if(A<1) call cq_abort("Error in matrix_pos: matrix number illegal: ",A)
    if(i>mat_p(A)%length.OR.i<1) call cq_abort("Matrix overrun in store_matrix_value_pos: ",i,mat_p(A)%length)
    !if(abs(val)>100.0_double) write(io_lun,*) 'Error: ',A,i,val
    ! In future we'll just have matrix_pos here
    mat_p(A)%matrix(i) = mat_p(A)%matrix(i) + val
!    call stop_timer(tmr_std_matrices)
    !if(abs(mat_p(A)%matrix(i))>100.0_double) write(io_lun,*) 'Error: ',A,i,val
  end subroutine store_matrix_value_pos

  subroutine store_matrix_block_pos(A,i,val,size)

    use datatypes
    use GenComms, ONLY: cq_abort
    use GenBlas, ONLY: axpy
    use numbers, ONLY: one

    implicit none

    ! Passed variables
    integer :: A, i, size
    real(double) :: val(size)

!    call start_timer(tmr_std_matrices)
    if(i+size-1>mat_p(A)%length) call cq_abort("Matrix overrun in store_matrix_value_pos: ",i+size-1,mat_p(A)%length)
    call axpy(size,one,val,1,mat_p(A)%matrix(i:),1)
!    call stop_timer(tmr_std_matrices)
  end subroutine store_matrix_block_pos

!  subroutine store_onsite_matrix_value(A,np,i,f1,f2,val)
!
!    use datatypes
!    use matrix_module, ONLY: matrix
!    use matrix_data, ONLY: mat
!
!    implicit none
!
!    ! Passed variables
!    integer :: A, np, i, f1, f2
!    real(double) :: val
!    ! Local variables
!    integer :: Ah, pos
!
!    Ah = matrix_index(A)
!    pos = mat(np,Ah)%onsite(i)
!    ! In future we'll just have matrix_pos here
!    mat_p(A)%matrix(f1,f2,pos) = val
!  end subroutine store_onsite_matrix_value

  integer function allocate_temp_matrix(range, trans_range, f1_type, f2_type)

    use matrix_data, ONLY: mat
    use GenComms, ONLY: cq_abort, myid
    use numbers, ONLY: zero
    use global_module, ONLY: iprint_mat, area_matrices
    use GenBlas, ONLY: scal
    use memory_module, ONLY: reg_alloc_mem, type_dbl

    implicit none

    ! Passed variables
    integer :: range, trans_range
    integer, optional :: f1_type, f2_type
    ! Local variables
    integer :: stat, i

    if(current_matrix==max_matrices) call cq_abort("Overrun number of matrices ",max_matrices)
    current_matrix = current_matrix + 1
    ! Identify the matrix
    if(range>mx_matrices) call cq_abort("Array error in alloc_temp_mat: ",range,mx_matrices)
    matrix_index(current_matrix) = range
    trans_index(current_matrix) = trans_range
    if(PRESENT(f1_type).AND.PRESENT(f2_type)) then
       mat_p(current_matrix)%sf1_type = f1_type
       mat_p(current_matrix)%sf2_type = f2_type
    else
       mat_p(current_matrix)%sf1_type = sf
       mat_p(current_matrix)%sf2_type = sf
    end if
    mat_p(current_matrix)%length = mat(1,range)%length
    call start_timer(tmr_std_allocation)
    allocate(mat_p(current_matrix)%matrix( mat_p(current_matrix)%length),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating matrix ",current_matrix,mat_p(current_matrix)%length)
    call stop_timer(tmr_std_allocation)
    !call scal(mat_p(current_matrix)%length,zero,mat_p(current_matrix)%matrix,1)
    do i=1,mat_p(current_matrix)%length
       mat_p(current_matrix)%matrix(i) = zero
    end do
    allocate_temp_matrix = current_matrix
    if(iprint_mat>4.AND.myid==0) write(io_lun,*) '(A)Current matrix is: ',current_matrix,mat_p(current_matrix)%length
    call reg_alloc_mem(area_matrices,mat_p(current_matrix)%length,type_dbl)
  end function allocate_temp_matrix

  subroutine free_temp_matrix(A)

    use GenComms, ONLY: cq_abort, myid
    use global_module, ONLY: iprint_mat, area_matrices
    use memory_module, ONLY: reg_dealloc_mem, type_dbl

    implicit none

    ! Passed variables
    integer :: A
    ! Local variables
    integer :: stat

    call reg_dealloc_mem(area_matrices,mat_p(current_matrix)%length,type_dbl)
    if(A/=current_matrix) call cq_abort("Out-of-order deallocation of matrices ",A,current_matrix)
    if(iprint_mat>4.AND.myid==0) write(io_lun,*) '(D)Current matrix is: ',current_matrix,mat_p(current_matrix)%length
    call start_timer(tmr_std_allocation)
    deallocate(mat_p(current_matrix)%matrix,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating matrix ",current_matrix,mat_p(current_matrix)%length)
    call stop_timer(tmr_std_allocation)
    mat_p(current_matrix)%length = 0
    matrix_index(current_matrix) = 0
    current_matrix = current_matrix - 1
  end subroutine free_temp_matrix

end module mult_module
