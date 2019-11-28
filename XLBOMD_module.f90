  !!****h* Conquest/XLBOMD_module
  !!  NAME
  !!   XLBOMD_module
  !!  PURPOSE
  !!   Performs the extended-Lagrangian Born-Oppenheimer MD 
  !!   adjusted to Conquest. Since Conquest deals with non-
  !!   orthogonal basis functions, it is necessary to modify
  !!   the original formulation.
  !!    - Original XL-BOMD: PRL 100, 123004 (2008)
  !!    - XL-BOMD with dissipation: J.Chem.Phys 130, 214109 (2009)
  !!  AUTHOR
  !!   Michiaki Arita
  !!  CREATION DATE
  !!   2013/12/03
  !!  MODIFICATION HISTORY
  !!   2017/05/11 dave
  !!    Various changes to allow L reuse and propagation with spin
  !!   2017/11/11 Tsuyoshi
  !!     io_module2 -> store_matrix
  !!     InfoX, InfoX1, ...., INfoXvel are defined in this module
  !!     InfoS is defined locally (in each subroutine).
  !!   2019/11/8  Tsuyoshi
  !!     grab_XXvelS, grab_Xhistories, dump_XL will(!) be moved to store_matrix,
  !!     since they will be called from the routines in store_matrix.
  !!     (this new module will be merged to store_matrix)
  !!   2019/11/15 tsuyoshi
  !!     1. matX, matXvel, matX_store are now moved to mult_module, and 
  !!     their names were changed to matXL, matXLvel, matXL_store.
  !!     2. the three routines (grab_XXvelS, ...) are now moved to store_matrix.
  !!  SOURCE
  !!
  module XLBOMD_module

    use datatypes
    use global_module, ONLY: iprint_MD,nspin
    use mult_module, ONLY: maxiter_Dissipation
    implicit none
    integer,allocatable :: matZ(:),matPot(:)
    !old integer,allocatable :: matX(:),matXvel(:),matZ(:),matPot(:)
    !old integer,allocatable :: matX_store(:,:)

    integer :: XLInitFreq
    real(double) :: kappa,kappa_diss,alpha
    real(double),allocatable :: c(:)

    character(80),private :: RCSid = "$Id$"
    logical, save :: allocated_XL = .false.

  contains

    ! ------------------------------------------------------------
    ! Subroutine immi_XL
    ! ------------------------------------------------------------

    !!****f* XLBOMD_module/immi_XL *
    !!  NAME
    !!   immi_XL 
    !!  USAGE
    !!   call immi_XL(parts,prim,gcs,myid)
    !!  PURPOSE
    !!   Allocates matrix X and its SCF solution Z (and its time derivative
    !!   when velocity verlet applies)
    !!   Initialises matrix multiplication type LS_T_L
    !!   Allocates matX_store to store and reorder a history of X-matrices
    !!   when dissipation applies
    !!  INPUTS
    !!   parts       : parts
    !!   prim        : bundle
    !!   gcs         : BCS_parts
    !!   integer myid: my process ID (myid+1 or inode)
    !!  AUTHOR
    !!   Michiaki Arita
    !!  CREATION DATE
    !!   2013/12/03 
    !!  MODIFICATION
    !!  SOURCE
    !!
    subroutine immi_XL(parts,prim,gcs,myid)
      ! Module usage
      use basic_types
      use numbers, ONLY: zero
      use global_module, ONLY: io_lun,integratorXL,flag_propagateL,flag_dissipation
      use GenComms, ONLY: cq_abort,inode,ionode
      use matrix_data, ONLY: mat,halo,LSrange,Lrange,Trange,Tmatind
      use mult_module, ONLY: mult,LS_trans,L_trans,LS_T_L,ltrans,allocate_temp_matrix
      use mult_module, ONLY: matXL, matXLvel, matXL_store,maxiter_Dissipation
      use mult_init_module, ONLY: mult_ini

      implicit none
      ! passed variables
      integer :: myid
      type(group_set),target   :: parts
      type(primary_set),target :: prim
      type(cover_set),target   :: gcs

      ! local variables
      integer :: ispin,i,stat_alloc
      integer :: range,trans,K
      logical,save :: allocated_tags1 = .false.
      logical,save :: allocated_tags2 = .false.
      ! db
      integer :: lun_db

      !if (inode.EQ.ionode) write (io_lun,*) "allocated_tags1 in XL: ", allocated_tags1
      !if (inode.EQ.ionode) write (io_lun,*) "allocated_tags2 in XL: ", allocated_tags2

      ! Allocate X, Z & Xvel
      if (.NOT.allocated_tags1) then
        allocate (matXL(nspin),matZ(nspin), STAT=stat_alloc)
        if (stat_alloc.NE.0) call cq_abort('Error allocating matXL and matZ: ', &
                                            nspin, stat_alloc)
        ! with velocity-Verlet
        if (integratorXL.EQ.'velocityVerlet') then
          allocate (matXLvel(nspin), STAT=stat_alloc)
          if (stat_alloc.NE.0) call cq_abort('Error allocating matXLvel: ', &
                                              nspin, stat_alloc)
        endif
        allocated_tags1 = .true.
      endif

      ! Choose the range
      if (.NOT. flag_propagateL) then
        range = LSrange
        trans = LS_trans
      else
        range = Lrange
        trans = L_trans
      endif
      do ispin = 1, nspin
        matXL(ispin)    = allocate_temp_matrix(range,trans)
        matZ(ispin)    = allocate_temp_matrix(range,trans)
        if (integratorXL.EQ.'velocityVerlet') &
           matXLvel(ispin) = allocate_temp_matrix(range,trans)
      enddo

      ! Initialise matrix multiplication
      if (.NOT. flag_propagateL) then
        mult(LS_T_L)%mult_type = 2
        mult(LS_T_L)%amat    => mat(1:prim%groups_on_node,Lrange)
        mult(LS_T_L)%bmat    => mat(1:prim%groups_on_node,Trange)
        mult(LS_T_L)%cmat    => mat(1:prim%groups_on_node,LSrange)
        mult(LS_T_L)%ahalo   => halo(Lrange)
        mult(LS_T_L)%chalo   => halo(LSrange)
        mult(LS_T_L)%ltrans  => ltrans(Lrange)
        mult(LS_T_L)%parts   => parts
        mult(LS_T_L)%prim    => prim
        mult(LS_T_L)%gcs     => gcs
        call mult_ini(mult(LS_T_L),Tmatind,myid-1,prim%n_prim,parts)
      endif

      ! When dissipation applies
      if (flag_dissipation) then
        K = maxiter_Dissipation
        if (.NOT.allocated_tags2) then
          allocate(matXL_store(K+1,nspin), STAT=stat_alloc)
          if (stat_alloc.NE.0) call cq_abort('Error allocating matXL_store: ', &
                                              stat_alloc)
          allocated_tags2 = .true.
        endif
        do ispin = 1, nspin
          do i = 1, K+1
            matXL_store(i,ispin) = allocate_temp_matrix(range,0)
          enddo
        enddo
      endif

      if (inode.EQ.ionode .AND. iprint_MD.GT.1) &
        write (io_lun,*) "Completed immi_XL()"

      return
    end subroutine immi_XL
    !!***

    ! ------------------------------------------------------------
    ! Subroutine  fmmi_XL
    ! ------------------------------------------------------------

    !!****f* XLBOMD_module/fmmi_XL *
    !!  NAME
    !!   fmmi_XL
    !!  USAGE
    !!   call fmmi_XL()
    !!  PURPOSE
    !!   Deallocates and nullifies all the matrix indexings for XL-BOMD
    !!  AUTHOR
    !!   Michiaki Arita
    !!  CREATION DATE
    !!   2013/12/03
    !!  MODIFICATION
    !!   2017/05/11 dave
    !!    Making consistent with immi_XL order with spin
    !!  SOURCE
    !!
    subroutine fmmi_XL
      ! Module usage
      use GenComms, ONLY: cq_abort,inode,ionode
      use global_module, ONLY: io_lun,integratorXL,flag_propagateL,flag_dissipation
      use matrix_module, ONLY: deallocate_comms_data
      use mult_module, ONLY: mult,LS_T_L,mult,free_temp_matrix
      use mult_module, ONLY: matXL_store, matXL, matXLvel, matXL_store

      implicit none
      ! local variables
      integer :: ispin,i,stat_alloc
      integer :: K

      !write (io_lun,*) "Got in deallocate_matXL", inode

      ! When dissipation applies
      if (flag_dissipation) then
        K = maxiter_Dissipation
        do ispin = nspin, 1, -1
          do i = K+1, 1, -1
            call free_temp_matrix(matXL_store(i,ispin))
          enddo
          !deallocate (matXL_store)
        enddo
      endif

      !if (allocated(matXLvel)) then
      do ispin = nspin,1,-1
         if (integratorXL.EQ.'velocityVerlet') &
              call free_temp_matrix(matXLvel(ispin))
         call free_temp_matrix(matZ(ispin))
         call free_temp_matrix(matXL(ispin))
      enddo
      !endif
      !if (allocated(matX) .AND. allocated(matZ)) then 
      !  do ispin = nspin,1,-1
      !  enddo
        ! Multiplication
      if (.NOT. flag_propagateL) then
         nullify(mult(LS_T_L)%amat,mult(LS_T_L)%bmat,mult(LS_T_L)%cmat,     &
              mult(LS_T_L)%ahalo,mult(LS_T_L)%chalo,mult(LS_T_L)%ltrans, &
              mult(LS_T_L)%bindex,mult(LS_T_L)%parts,mult(LS_T_L)%prim,  &
              mult(LS_T_L)%gcs)
         ! Comms
         call deallocate_comms_data(mult(LS_T_L)%comms)
      endif
      !endif

      if (inode.EQ.ionode .AND. iprint_MD.GT.1) &
        write (io_lun,*) "Completed fmmi_XL()"

      return
    end subroutine fmmi_XL
    !!***

    ! ------------------------------------------------------------
    ! Subroutine Do_XLBOMD
    ! ------------------------------------------------------------

    !!****f* XLBOMD_module/Do_XLBOMD *
    !!  NAME
    !!   Do_XLBOMD
    !!  USAGE
    !!   call Do_XLBOMD(MDiter,dt)
    !!  PURPOSE
    !!   Reads and reconstructs matrices used in XL-BOMD then
    !!   propagates X-matrix
    !!  INPUTS
    !!   integer MDiter: current MD step
    !!   real(double) dt: MD time step
    !!  AUTHOR
    !!   Michiaki Arita
    !!  CREATION DATE
    !!   2013/12/03
    !!  MODIFICATION
    !!  SOURCE
    subroutine Do_XLBOMD(MDiter,dt)
      ! Module usage
      use global_module, ONLY: flag_propagateX,flag_propagateL
      use matrix_data, ONLY: LSrange,Lrange
      use mult_module, ONLY: LS_trans,L_trans,matL

      implicit none
      ! passed variables
      integer,intent(in) :: MDiter
      real(double),intent(in) :: dt

      if (flag_propagateX) then
        call CommRebuild_matXL(MDiter,LSrange,LS_trans)
        call propagate_matXL(MDiter,dt,matZ,LSrange)
      elseif (flag_propagateL) then
        call CommRebuild_matXL(MDiter,Lrange,L_trans)
        call propagate_matXL(MDiter,dt,matL,Lrange)
      endif

      return
    end subroutine Do_XLBOMD
    !!***

    ! ------------------------------------------------------------
    ! Subroutine get_initialL_XL
    ! ------------------------------------------------------------

    !!****f* XLBOMD_module/get_initialL_XL *
    !!  NAME
    !!   get_initialL_XL
    !!  USAGE
    !!   call get_initialL_XL()
    !!  PURPOSE
    !!   Performs;
    !!      1. L^{init}(t+dt) = X(t+dt)*T(t+dt) for modified scheme
    !!      2. L^{init}(t+dt) = X(t+dt)         for original scheme
    !!  AUTHOR
    !!   Michiaki Arita
    !!  CREATION DATE
    !!   2013/12/03
    !!  MODIFICATION
    !!   2018/11/13 17:30 nakata
    !!    Changed matT and matTtran to be spin_SF dependent
    !!  SOURCE
    !!
    subroutine get_initialL_XL()
      ! Module usage
      use numbers
      use global_module, ONLY: nspin,io_lun,flag_propagateX,flag_propagateL, &
                               nspin_SF, flag_SpinDependentSF
      use GenComms, ONLY: inode,ionode
      use mult_module, ONLY: mult,matrix_product,matL,matT,matTtran,LS_T_L, &
                             symmetrise_L,mat_p,matrix_sum,matrix_transpose, matXL

      implicit none
      ! local variables
      integer :: ispin, ispin_SF

      if (inode.EQ.ionode .AND. iprint_MD.GT.1) &
        write (io_lun,*) "Make an initial guess on L-matrix"

      if (flag_propagateX) then
        call matrix_transpose(matT(1),matTtran(1))
        if (nspin_SF==2) call matrix_transpose(matT(2),matTtran(2))
        ispin_SF = 1
        do ispin = 1, nspin
          if (flag_SpinDependentSF) ispin_SF = ispin
          !WRONG call matrix_product(matXL(ispin),matT(ispin_SF),matL(ispin),mult(LS_T_L))
          call matrix_product(matXL(ispin),matTtran(ispin_SF),matL(ispin),mult(LS_T_L))
        enddo
      elseif (flag_propagateL) then
        !ORI mat_p(matL(1))%matrix = mat_p(matXL(1))%matrix
        do ispin = 1, nspin
          call matrix_sum(zero,matL(ispin),one,matXL(ispin))
        enddo
      endif
      call symmetrise_L()

      return
    end subroutine get_initialL_XL
    !!***

    ! ------------------------------------------------------------
    ! Subroutine CommRebuild_matXL
    ! ------------------------------------------------------------

    !!****f* XLBOMD_module/CommRebuild_matXL *
    !!  NAME
    !!   CommRebuild_matXL
    !!  USAGE
    !!   call CommRebuild_matXL(MDiter,range,trans)
    !!  PURPOSE
    !!   Calculates Z(t) = L^{opt}(t)*S(t).
    !!   First grabs L,S,X and Xvel matrices and reconstructs them.
    !!   Then performs matrix multiplication.
    !!   This subroutine must be called before calculating S at
    !!   t=t+dt, or calling update_H. Also handles initialisation
    !!   and storing & reordering X-matrices
    !!  INPUTS
    !!   integer MDiter: current MD step
    !!   integer range : for modified ver. LSrange, otherwise Lrange
    !!   integer trans : for modified ver. LS_trans, otherwise L_trans
    !!  AUTHOR
    !!   Michiaki Arita
    !!  CREATION DATE
    !!   2013/12/03
    !!  MODIFICATION
    !!   2018/11/13 17:30 nakata
    !!    Changed matS to be spin_SF dependent
    !!  SOURCE
    !!
    subroutine CommRebuild_matXL(MDiter,range,trans)
      ! Module usage
      use global_module, ONLY: nspin,flag_dissipation,integratorXL,flag_propagateX, &
                               flag_SpinDependentSF
      use GenComms, ONLY: cq_abort,inode,my_barrier
      use matrix_data, ONLY: LSrange,Srange
      use mult_module, ONLY: LS_trans,S_trans,matS,matL,matrix_product,mult, &
                             L_S_LS
      ! db
      use global_module, ONLY: io_lun
      use GenComms, ONLY: myid
 !2019/Nov/13
      use numbers, ONLY: half, one, very_small
      use global_module, ONLY: io_lun,ni_in_cell,atom_coord,atom_coord_diff,rcellx,rcelly,rcellz
      use GenComms, ONLY: ionode
      use store_matrix,ONLY: matrix_store_global, grab_InfoMatGlobal, set_atom_coord_diff

      implicit none
      ! passed variables
      integer,intent(in) :: MDiter,range,trans
      ! local variables
      integer :: nfile,ispin,ispin_SF,symm
 !2019/Nov/13
      type(matrix_store_global) :: InfoGlob   

 !2019/Nov/13 InfoGlob
      call grab_InfoMatGlobal(InfoGlob,index=0)
      call set_atom_coord_diff(InfoGlob)

      ! Fetches & reconstructs X-matrix, Xvel-matrix & S-matrix
      call grab_XXvelS(range,trans,InfoGlob)
      ! Gets Z(t) = L^{opt}(t) * S(t)
      if (flag_propagateX) then
        ispin_SF = 1
        do ispin = 1, nspin
          if (flag_SpinDependentSF) ispin_SF = ispin
          call matrix_product(matL(ispin),matS(ispin_SF),matZ(ispin),mult(L_S_LS))
        enddo
        if (myid.EQ.0 .AND. iprint_MD.GT.1) &
          write (io_lun,*) "Got Z-matrix" 
      endif

      ! Initialises or resets X-matrix
      call init_matXL()

      if (flag_dissipation) then
        call grab_Xhistories(range,trans,InfoGlob)
        call my_barrier()
        call reorder_Xhistories(MDiter)
        if (myid.EQ.0 .AND. iprint_MD.GT.1) &
          write (io_lun,*) "Completed reorder_Xhistories"
      endif

      return
      contains

      subroutine init_matXL()
        ! Module usage
        use numbers, ONLY: zero,one
        use global_module, ONLY: integratorXL,flag_propagateX,flag_SpinDependentSF
        use GenComms, ONLY: cq_abort
        use mult_module, ONLY: mult,matrix_product,matrix_sum,matS,matL,L_S_LS, &
                               mat_p,matrix_scale, matXL, matXLvel

        implicit none
        ! local variables
        integer,save :: iter = 0
        integer :: ispin,ispin_SF,n
        logical :: flag

        flag = .false.
        iter = iter + 1
        if (XLInitFreq.NE.0) then
          n = mod(iter,XLInitFreq)
        else
          n = 1
        endif
        if (iter.EQ.1 .OR. n.EQ.0) flag = .true.
        if (flag) then
          if (myid.EQ.0 .AND. iprint_MD.GT.1) &
            write (io_lun,*) "Initialising X-matrix", iter
          if (flag_propagateX) then
            ispin_SF = 1
            do ispin = 1, nspin
              if (flag_SpinDependentSF) ispin_SF = ispin
              call matrix_product(matL(ispin),matS(ispin_SF),matXL(ispin),mult(L_S_LS))
              call matrix_sum(zero,matZ(ispin),one,matXL(ispin))
!% michi - comments 14/11/2013
!%            I think the above two calls are redundant. It should
!%            call matrix_sum(zero,matXL(ispin),one,matZ(ispin))
!%            because matZ is already obtained before this subroutine gets called.
!% michi - comments 14/11/2013
              if (integratorXL.EQ.'velocityVerlet') &
                call matrix_scale(zero,matXLvel(ispin))
            enddo
          else
            do ispin = 1, nspin
              call matrix_sum(zero,matXL(ispin),one,matL(ispin))
              if (integratorXL.EQ.'velocityVerlet') &
                call matrix_scale(zero,matXLvel(ispin))
            enddo
          endif
        endif

        return
      end subroutine init_matXL

    end subroutine CommRebuild_matXL
    !!***

    ! ------------------------------------------------------------
    ! Subroutine grab_XXvelS
    ! ------------------------------------------------------------

    !!****f* XLBOMD_module/grab_XXvelS *
    !!  NAME
    !!   grab_XXvelS
    !!  USAGE
    !!   call grab_XXvelS(range,trans)
    !!  PURPOSE
    !!   Reads and reconstructs X-matrix, Xvel-matrix and S-matrix
    !!  INPUTS
    !!   integer range: for modified ver. LSrange, otherwise Lrange
    !!   integer trans: for modified ver. LS_trans, otherwise L_trans
    !!  AUTHOR
    !!   Michiaki Arita
    !!  CREATION DATE
    !!   2013/12/03 
    !!  MODIFICATION
    !!   2017/05/11 dave
    !!    Adding read option for spin polarisation
    !!   2018/11/13 17:30 nakata
    !!    Changed matS to be spin_SF dependent
    !!   2018/11/15 15:45 nakata
    !!    Bug fix: InfoXvel should be used for Xvel_2
    !!  SOURCE
    subroutine grab_XXvelS(range,trans,InfoGlob)
      ! Module usage
      use global_module, ONLY: integratorXL,flag_propagateX, nspin, nspin_SF
      use GenComms, ONLY: inode, ionode
      use matrix_data, ONLY: Srange
      use mult_module, ONLY: matS,S_trans,matXL,matXLvel
      use store_matrix, ONLY: grab_matrix2, InfoMatrixFile, matrix_store_global
      use UpdateInfo, ONLY: Matrix_CommRebuild

      implicit none
!      ! passed variables
      integer,intent(in) :: range,trans
      ! local variables
      integer :: nfile,symm
      type(InfoMatrixFile), pointer :: InfoMat(:)
 !2019/Nov/08
      type(matrix_store_global),intent(in) :: InfoGlob   

      ! Fetches & reconstructs X-matrix
      call grab_matrix2('X',inode,nfile,InfoMat,InfoGlob,index=0,n_matrix=nspin)
      !call Matrix_CommRebuild(InfoGlob,InfoMat,range,trans,matXL,nfile,symm,n_matrix=nspin)
      call Matrix_CommRebuild(InfoGlob,InfoMat,range,trans,matXL,nfile,n_matrix=nspin)
      ! Fetches & reconstructs Xvel-matrix
      if (integratorXL.EQ.'velocityVerlet') then
        call grab_matrix2('Xvel',inode,nfile,InfoMat,InfoGlob,index=0,n_matrix=nspin)
        !call Matrix_CommRebuild(InfoGlob,InfoMat,range,trans,matXLvel,symm,nfile,n_matrix=nspin)
        call Matrix_CommRebuild(InfoGlob,InfoMat,range,trans,matXLvel,nfile,n_matrix=nspin)
      endif
      ! Fetches & reconstructs S-matrix
      if (flag_propagateX) then
        call grab_matrix2('S',inode,nfile,InfoMat,InfoGlob,index=0,n_matrix=nspin_SF)
        call Matrix_CommRebuild(InfoGlob,InfoMat,Srange,S_trans,matS,nfile,symm,n_matrix=nspin_SF)
      endif

      return
    end subroutine grab_XXvelS
    !!***
!
    ! ------------------------------------------------------------
    ! Subroutine grab_Xhistories
    ! ------------------------------------------------------------

    !!****f* XLBOMD_module/grab_Xhistories *
    !!  NAME
    !!   grab_Xhistories
    !!  USAGE
    !!   call grab_Xhistories(range,trans)
    !!  PURPOSE
    !!   Reads and reconstructs a history of X-matrices
    !!  INPUTS
    !!   integer range: for modified ver. LSrange, otherwise Lrange
    !!   integer trans: for modified ver. LS_trans, otherwise L_trans
    !!  AUTHOR
    !!   Michiaki Arita
    !!  CREATION DATE
    !!   2013/12/03 
    !!  MODIFICATION
    !!   2017/05/11 dave
    !!    Adding read option for spin polarisation
    !!   2019/05/24 tsuyoshi
    !!    Changing the filenames and tidying up the code
    !!   2019/11/15 tsuyoshi
    !!   
    !!  SOURCE
    subroutine grab_Xhistories(range,trans,InfoGlob)
      ! Module usage
      use global_module, ONLY: flag_propagateL, nspin, io_lun
      use GenComms, ONLY: inode, ionode 
      use matrix_data, ONLY: LSrange,Lrange
      use mult_module, ONLY: LS_trans,L_trans, matXL_store, maxiter_Dissipation
      use store_matrix, ONLY: grab_matrix2, InfoMatrixFile, matrix_store_global
      use UpdateInfo, ONLY: Matrix_CommRebuild
      !db
      use global_module, ONLY: io_lun
      use GenComms, ONLY: myid

      implicit none
      ! passed variables
      integer :: range,trans
      !2019/Nov/13
      type(matrix_store_global),intent(in) :: InfoGlob   
      !integer :: matX_store(maxiter_Dissipation,nspin)
      ! local variables
      integer :: maxiters,nfile,istep
      type(InfoMatrixFile), pointer :: InfoMat(:)

      maxiters = maxiter_Dissipation
      if(maxiters < 3) maxiters=3  ! even without dissipation, we need matX_store(:,1:4)
      if(maxiters > 9) then
       if(inode == ionode) write(io_lun,*)  &
        &'WARNING: maxiter_Dissipation should be smaller than 10 : ', maxiter_Dissipation
       maxiters = 9 
      endif

      ! ----- 2019/Nov/13: (comment by TM)   -----
      ! Here, we assume InfoGlob is read before calling this subroutine.
      ! We use same InfoGlob for all index, meaning that the matrices stored in Xmatrix.i**,p**** are
      ! already rebuilt using the pair information shown in InfoGlob.i00, the one at the last step)
      !

      ! Grab X-matrix files
      do istep = 1, maxiters+1
       call grab_matrix2('X',inode,nfile,InfoMat,InfoGlob,index=istep,n_matrix=nspin)
       call Matrix_CommRebuild(InfoGlob,InfoMat,range,trans,matXL_store(istep,:),nfile,n_matrix=nspin)
      enddo

      ! ----- 2019/Nov/13: (comment by TM)   -----
      ! 1. Use same index in InfoGlobal.iXX and Xmatrix2.iXX.pYYYYYY
      ! 2. Use InfoGlobal.i00 for all Xmatrix2.iXX.pYYYYYY   
      ! 
      !  Present scheme uses Option 2, which should be more efficient.
      ! But, in this case, it is more efficient if we call Matrix_CommRebuild once, 
      ! with  n_matrix = (maxiters+1)*nspin.   
      !  (since it is expensive to make the transformation table (i,j) <-> (i_old,j_old))
      !
      ! When we use Memory for InfoMat(:) and InfoGlob, instead of using IOs, we should prepare
      !  set_prev_matrix for grab_matrix2, and we should be able to rebuild matX_store once,
      ! by setting n_matix = (maxiters+1)*nspin
      !

      return
    end subroutine grab_Xhistories
    !!***
!
!    ! ------------------------------------------------------------
!    ! Subroutine dump_XL
!    ! ------------------------------------------------------------
!
!    !!****f* XLBOMD_module/dump_XL *
!    !!  NAME
!    !!   dump_XL
!    !!  USAGE
!    !!   call dump_XL()
!    !!  PURPOSE
!    !!   Dumps all matrix files related to XL-BOMD
!    !!  AUTHOR
!    !!   Michiaki Arita
!    !!  CREATION DATE
!    !!   2013/12/03 
!    !!  MODIFICATION
!    !!   2017/05/11 dave
!    !!    Adding write option for spin polarisation
!    !!   2019/05/24 tsuyoshi
!    !!    Change the filenames and tidying up the code
!    !!  SOURCE
!    subroutine dump_XL()
!      ! Module usage
!      use global_module, ONLY: flag_propagateL, nspin,io_lun
!      use GenComms, ONLY: cq_abort,inode,ionode
!      use matrix_data, ONLY: LSrange,Lrange
!      use store_matrix, ONLY: dump_matrix2
!
!      implicit none
!      ! local variables
!      integer :: nfile,maxiters,range,istep
!
!      if (.NOT. flag_propagateL) then
!        range = LSrange
!      else
!        range = Lrange
!      endif
!      maxiters = maxiter_Dissipation
!      if(maxiters < 4) maxiters = 3
!      if(maxiters > 9) then
!       if(inode == ionode) write(io_lun,*)  &
!        &'WARNING: maxiter_Dissipation should be smaller than 10 : ', maxiter_Dissipation
!       maxiters = 9 
!      endif
!
!      ! Dump X-matrix files
!      do istep = 1, maxiters+1
!        !call dump_matrix2('X',matX_store(istep,1),range,index=istep)
!        call dump_matrix2('X',matX_store(istep,:),range,n_matrix=nspin,index=istep)
!      enddo
!      !if(nspin==2) then
!      ! do istep = 1, maxiters+1
!      !  call dump_matrix2('X2',matX_store(istep,2),range,index=istep)
!      ! enddo
!      !endif
!
!      return
!    end subroutine dump_XL
!    !!***
!
    ! ------------------------------------------------------------
    ! Subroutine reorder_Xhistories
    ! ------------------------------------------------------------

    !!****f* XLBOMD_module/reorder_XLhistories *
    !!  NAME
    !!   reorder_Xhistories
    !!  USAGE
    !!   call reorder_Xhistories(MDiter)
    !!  PURPOSE
    !!   For MDiter .GT. K+1, reorders X-matrix files, otherwise stores
    !!   X-matrix in X?matrix2.???. The oldest X-matrix is written in 
    !!   X1matrix2.???.
    !!  INPUT
    !!   integer MDiter: current MD step
    !!  AUTHOR
    !!   Michiaki Arita
    !!  CREATION DATE
    !!   2013/12/03 
    !!  MODIFICATION
    !!  SOURCE
    !!
    subroutine reorder_Xhistories(MDiter)
      ! Module usage
      use numbers
      use mult_module, ONLY: matrix_sum, matXL_store, matXL
      ! db
      use global_module, ONLY: io_lun
      use GenComms, ONLY: myid

      implicit none
      ! passed variables
      integer,intent(in) :: MDiter
      ! local variables
      integer :: i,ispin
      integer :: K

      K = maxiter_Dissipation
      if (MDiter.GT.K+1) then
        if (myid.EQ.0 .AND. iprint_MD.GT.2) &
          write (io_lun,*) "Reordering X-matrix files",MDiter
        do ispin = 1,nspin
          do i = 1, K
            call matrix_sum(zero,matXL_store(i,ispin),one,matXL_store(i+1,ispin))
          enddo
          call matrix_sum(zero,matXL_store(K+1,ispin),one,matXL(ispin))
        enddo
      elseif (MDiter.LE.K+1) then
        if (myid.EQ.0 .AND. iprint_MD.GT.2) &
          write (io_lun,*) "Storing X-matrix",MDiter
        do ispin = 1, nspin
          call matrix_sum(zero,matXL_store(MDiter,ispin),one,matXL(ispin))
        enddo
      endif

      return
    end subroutine reorder_Xhistories
    !!***

    ! ------------------------------------------------------------
    ! Subroutine calc_dissipative_force
    ! ------------------------------------------------------------

    !!****f* XLBOMD_module/calc_dissipative_force *
    !!  NAME
    !!   calc_dissipative_force
    !!  USAGE
    !!   call calc_dissipative_force(MDiter,dissipation)
    !!  PURPOSE
    !!   Calculates a dissipative force to eliminate uknown
    !!   inherent numerical noise
    !!     dissipation = alpha * sum_{k=1}^{K+1} c_k * X(t-(k-1)*dt)
    !!  INPUTS
    !!   integer MDiter: current MD step
    !!   real(double),dimension(K+1,nspin) dissipation: dissipative force
    !!  AUTHOR
    !!   Michiaki Arita
    !!  CREATION DATE
    !!   2013/12/03 
    !!  MODIFICATION
    !!  SOURCE
    !!
    subroutine calc_dissipative_force(MDiter,dissipation)
      ! Module usage
      use numbers, ONLY: zero
      use mult_module, ONLY: mat_p, matXL, matXL_store
      !db
      use global_module, ONLY: numprocs,io_lun
      use GenComms, ONLY: inode,myid
      use input_module, ONLY: io_assign,io_close
      use io_module, ONLY: get_file_name

      implicit none
      ! passed variables
      integer,intent(in) :: MDiter
      real(double) :: dissipation(:,:)
      ! local variables
      integer :: i,ispin,len,K
      !db
      integer :: lun_db
      character(20) :: file_name


      !db if (myid.EQ.0) write (io_lun,'(a,1x,i8,1x,a)') &
      !db                      "Calculate dissipative force at", MDiter, "MD iter"

      dissipation = zero
      K = maxiter_Dissipation
      do ispin = 1, nspin
        len=mat_p(matXL(ispin))%length
        do i = 1, K+1
            dissipation(1:len,ispin) = dissipation(1:len,ispin) + &
                                     c(i)*mat_p(matXL_store(i,ispin))%matrix(1:len)
        enddo
      enddo
      dissipation = dissipation * alpha

      return
    end subroutine calc_dissipative_force
    !!***

    ! ------------------------------------------------------------
    ! Subroutine Ready_XLBOMD
    ! ------------------------------------------------------------

    !!****f* XLBOMD_module/Ready_XLBOMD *
    !!  NAME
    !!   Redy_XLBOMD
    !!  USAGE
    !!   call Ready_XLBOMD
    !!  PURPOSE
    !!   Gets optimal parameters for a dissipative force
    !!   Refer to J.Chem.Phys 130, 214109 (2009)
    !!  AUTHOR
    !!   Michiaki Arita
    !!  CREATION DATE
    !!   2013/12/03 
    !!  MODIFICATION
    !!  SOURCE
    !!
    subroutine Ready_XLBOMD()
      ! Module usage
      use global_module, ONLY: io_lun
      use GenComms, ONLY: cq_abort,myid

      implicit none
      ! local variables
      integer :: i,stat_alloc,K

      K = maxiter_Dissipation
      allocate (c(K+1), STAT=stat_alloc)
      if (stat_alloc.NE.0) call cq_abort('Error allocating coefficcient c: ', &
                                         K+1)

      !%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
      !%    NOTE:					     %!
      !%      matX_store(1:K+1,ispin) & c(1:K+1)	     %!
      !%        - K <-- newer , older --> 1                  %!
      !%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!

      select case (K)
      case (3)
        kappa_diss = 1.69_double
        alpha = 0.15_double
        c(1)  = -1.0_double
        c(2)  = 0.0_double
        c(3)  = 3.0_double
        c(4)  = -2.0_double
      case (4)
        kappa_diss = 1.75_double
        alpha = 0.057_double
        c(1)  = 1.0_double
        c(2)  = -2.0_double
        c(3)  = -2.0_double
        c(4)  = 6.0_double
        c(5)  = -3.0_double
      case (5)
        kappa_diss = 1.82_double
        alpha = 0.018_double
        c(1)  = -1.0_double
        c(2)  = 4.0_double
        c(3)  = -3.0_double
        c(4)  = -8.0_double
        c(5)  = 14.0_double
        c(6)  = -6.0_double
      case (6)
        kappa_diss = 1.84_double
        alpha = 0.0055_double
        c(1)  = 1.0_double
        c(2)  = -6.0_double
        c(3)  = 12.0_double
        c(4)  = -2.0_double
        c(5)  = -27.0_double
        c(6)  = 36.0_double
        c(7)  = -14.0_double
      case (7)
        kappa_diss = 1.86_double
        alpha = 0.0016_double
        c(1)  = -1.0_double
        c(2)  = 8.0_double
        c(3)  = -25.0_double
        c(4)  = 32.0_double
        c(5)  = 11.0_double
        c(6)  = -88.0_double
        c(7)  = 99.0_double
        c(8)  = -36.0_double
      case (8)
        kappa_diss = 1.88_double
        alpha = 0.00044_double
        c(1)  = 1.0_double
        c(2)  = -10.0_double
        c(3)  = 42.0_double
        c(4)  = -90.0_double
        c(5)  = 78.0_double
        c(6)  = 78.0_double
        c(7)  = -286.0_double
        c(8)  = 286.0_double
        c(9)  = -99.0_double
      case (9)
        kappa_diss = 1.89_double
        alpha = 0.00012_double
        c(1)  = -1.0_double
        c(2)  = 12.0_double
        c(3)  = -63.0_double
        c(4)  = 184.0_double
        c(5)  = -300.0_double
        c(6)  = 168.0_double
        c(7)  = 364.0_double
        c(8)  = -936.0_double
        c(9)  = 858.0_double
        c(10) = -286.0_double
      end select

      ! Corrects kappa if necessary
      if (kappa.NE.kappa_diss) then
        if (myid.EQ.0) write (io_lun,'(a,f8.5)') &
                       "WARNING: kappa is set to ",kappa_diss
      endif

      return
    end subroutine Ready_XLBOMD
    !!***

    ! ------------------------------------------------------------
    ! Subroutine propagate_matXL
    ! ------------------------------------------------------------

    !!****f* XLBOMD_module/propagate_matXL *
    !!  NAME
    !!   propagate_matXL
    !!  USAGE
    !!   call propagate_matXL(MDiter,dt,matA,Arange)
    !!  PURPOSE
    !!   Propagates X-matrix either with velocity Verlet or Verlet
    !!  INPUTS
    !!   integer MDiter: current MD step
    !!   real(double) dt: time step
    !!   integer matA: A-matrix tag
    !!   integer Arange: A-matrix range
    !!  AUTHOR
    !!   Michiaki Arita
    !!  CREATION DATE
    !!   2013/12/03 
    !!  MODIFICATION
    !!  SOURCE
    !!
    subroutine propagate_matXL(MDiter,dt,matA,Arange)
      ! Module usage
      use global_module, ONLY: integratorXL
      use GenComms, ONLY: cq_abort

      implicit none
      ! passed variables
      integer,intent(in) :: MDiter,matA(nspin),Arange
      real(double),intent(in) :: dt

      if (integratorXL.EQ.'Verlet') then
        call Verlet_matXL(MDiter,matA,Arange)
      elseif (integratorXL.EQ.'velocityVerlet') then
        call vVerlet_matXL(MDiter,dt,matA,Arange)
      else
        call cq_abort('Error: choose either Verlet or velocityVerlet')
      endif

      return
    end subroutine propagate_matXL
    !!***

    ! ------------------------------------------------------------
    ! Subroutine Verlet_matXL
    ! ------------------------------------------------------------

    !!****f* XLBOMD_module/Verlet_matXL *
    !!  NAME
    !!   Verlet_matXL
    !!  USAGE
    !!   call Verlet_matXL(MDiter,matA,Arange)
    !!  PURPOSE
    !!   Performs the Verlet algorithm for propagation of X-matrix
    !!     X(t+dt) = 2X(t) - X(t-dt) + (kappa)*[Z(t)-X(t)]
    !!  INPUTS
    !!   integer MDiter: current MD step
    !!   integer matA: A-matrix tag
    !!   integer Arange: A-matrix range
    !!  AUTHOR
    !!   Michiaki Arita
    !!  CREATION DATE
    !!   2013/12/03 
    !!  MODIFICATION
    !!   2018/11/13 17:30 nakata
    !!    Removed matS, matL, L_S_LS and mult which was passed but not used
    !!  SOURCE
    !!
    subroutine Verlet_matXL(MDiter,matA,Arange)
      ! Module usage
      use numbers
      use global_module, ONLY: flag_dissipation
      use GenComms, ONLY: inode,cq_abort
      use matrix_data, ONLY: LSrange
      use mult_module, ONLY: mat_p,matrix_product, &
                             matrix_sum,allocate_temp_matrix,free_temp_matrix, &
                             matrix_product_trace, &
                             matXL,matXL_store
      ! db
      use global_module, ONLY: io_lun
      use GenComms, ONLY: myid

      implicit none
      ! passed variables
      integer,intent(in) :: MDiter,matA(nspin),Arange
      ! local variables
      integer :: len,ispin,stat,j,K
      real(double) :: kappa_wk,acc,Pot_el,dr,dr_dt
      real(double),allocatable :: Fdiss(:,:)
      logical,save :: allocated_tag = .false.

      !% NOTE: In Verlet, X(dt) = X(0) + dt*Xvel(0) + half*kappa*(Z(0)-X(0))
      !%                        = X(0)
      !%       , so no need to integrate at MDiter = 1
      if (MDiter.EQ.1) return

      ! Allocate matPot
      if (.NOT.allocated_tag) then
        allocate (matPot(nspin), STAT=stat)
        if (stat.NE.0) call cq_abort('Error allocating matPot: ', nspin)
        allocated_tag = .true.
      endif
      do ispin = 1, nspin
        matPot(ispin) = allocate_temp_matrix(Arange,0)
      enddo

      if (MDiter.GT.1) then
        if (flag_dissipation) then
          K = maxiter_Dissipation
          allocate (Fdiss(mat_p(matXL(1))%length,nspin), STAT=stat)
          if (stat.NE.0) call cq_abort('Error allocating dissipation: ')
          Fdiss = zero
          if (MDiter.GT.K) call calc_dissipative_force(MDiter,Fdiss)
          ! Decide 'j'
          if (MDiter.LE.K+1) then
            j = MDiter-1
          else
            j = K
          endif
          kappa_wk = kappa_diss
          do ispin = 1, nspin
            !%ORI do len = 1, mat_p(matXL(1))%length
            !%ORI   ! electronic harmonic potential
            !%ORI   acc = mat_p(matZ(ispin))%matrix(len) - mat_p(matXL(ispin))%matrix(len)
            !%ORI   mat_p(matPot(ispin))%matrix(len) = acc
            !%ORI   mat_p(matXL(ispin))%matrix(len) = two*mat_p(matXL(ispin))%matrix(len)     - &
            !%ORI                                    mat_p(matXL_store(j,ispin))%matrix(len) + &
            !%ORI                                    kappa_wk*acc + Fdiss(len,ispin)
            !%ORI enddo
            do len = 1, mat_p(matXL(ispin))%length
              ! electronic harmonic potential
              acc    = mat_p(matA(ispin))%matrix(len) - mat_p(matXL(ispin))%matrix(len)
              mat_p(matPot(ispin))%matrix(len) = acc
              dr     = mat_p(matXL(ispin))%matrix(len) - mat_p(matXL_store(j,ispin))%matrix(len)
              dr_dt  = dr + kappa_wk*acc + Fdiss(len,ispin)
              mat_p(matXL(ispin))%matrix(len) = mat_p(matXL(ispin))%matrix(len) + dr_dt
            enddo
          enddo
          ! Deallocation
          deallocate (Fdiss, STAT=stat)
          if (stat.NE.0) call cq_abort('Error deallocating dissipation: ')
        else
          call cq_abort('Error: without dissipation, Verlet is unavailable !')
        endif
      endif
      ! Calculate harmonic terms
      Pot_el = half*matrix_product_trace(matPot(1),matPot(1))
      if (myid.EQ.0) write (io_lun,'(a,2x,i8,2x,a,f25.15)') &
                           "     *** MD iter",MDiter,'Tr[(Z-X)^2]:',Pot_el

      ! Free matPot
      do ispin = nspin,1,-1
        call free_temp_matrix(matPot(ispin))
      enddo

      if (myid.EQ.0 .AND. iprint_MD.GT.1) &
        write (io_lun,*) "X-matrix propagated via Verlet"

      return
    end subroutine Verlet_matXL
    !!***

    ! ------------------------------------------------------------
    ! Subroutine vVerlet_matXL
    ! ------------------------------------------------------------

    !!****f* XLBOMD_module/vVerlet_matXL *
    !!  NAME
    !!   vVerlet_matXL
    !!  USAGE
    !!   call vVerlet_matXL(MDiter,dt,matA,Arange)
    !!  PURPOSE
    !!   Performs the velocity Verlet algorithm for propagating X-matrix
    !!   The defifnition of velocities is a little tricky:
    !!     Xvel_init = Xvel_init(-dt/2)
    !!
    !!      1. Xvel(t)      = Xvel(t-dt/2) + (kappa/(2*dt))*[Z(t)-X(t)]
    !!      2. X(t+dt)      = X(t) + dt*Xvel(t) + (kappa/2)*[Z(t)-X(t)]
    !!      3. Xvel(t+dt/2) = Xvel(t) + (kappa/(2*dt))*[Z(t)-X(t)]
    !!
    !!   Accelaration becomes as follows with dissipation
    !!     ¥ddot{X}(t) = omega^2 * [Z(t) - X(t)] + (alpha/dt^2) * ¥sum_{k=1}^{K+1} X(t+(k-1)*dt)
    !!                 = omega^2 * [Z(t) - X(t)] + (1/dt^2) * dissipation
    !!
    !!  INPUTS
    !!   integer MDiter: current MD step
    !!   real(double) dt: time step
    !!   integer matA: matrix tag
    !!   integer Arange: A-matrix range
    !!  AUTHOR
    !!   Michiaki Arita
    !!  CREATION DATE
    !!   2013/12/03 
    !!  MODIFICATION
    !!  SOURCE
    !!
    subroutine vVerlet_matXL(MDiter,dt,matA,Arange)
      ! Module usage
      use numbers
      use global_module, ONLY: io_lun,flag_dissipation
      use GenComms, ONLY: cq_abort,myid
      use matrix_data, ONLY: LSrange
      use mult_module, ONLY: mat_p,matrix_product_trace,allocate_temp_matrix, &
                             free_temp_matrix, matXL, matXLvel

      implicit none
      ! passed variables
      integer     ,intent(in) :: MDiter,matA(nspin),Arange
      real(double),intent(in) :: dt
      ! local variables
      integer :: ispin,LSmatsize,len,stat_alloc
      real(double) :: acc,KE_el,kappa_wk,acc_diss,Pot_el
      real(double),allocatable :: Fdiss(:,:)
      logical,save :: allocated_tag = .false.
      !db
      integer :: i,n

      ! Allocate matPot
      if (.NOT.allocated_tag) then
        allocate (matPot(nspin), STAT=stat_alloc)
        if (stat_alloc.NE.0) call cq_abort('Error allocating matPot:')
        allocated_tag = .true.
      endif
      do ispin = 1, nspin
        matPot(ispin) =  allocate_temp_matrix(Arange,0)
      enddo

      ! Evolve matXL
      if (flag_dissipation) then

        !db
        !if (myid.EQ.0) then
        !  write (io_lun,'(a,f10.5)') "kappa:  ", kappa_diss
        !  write (io_lun,'(a,f10.5)') "alpha:  ", alpha
        !  write (io_lun,'(a)')       "c(1:K): "
        !  do i = 1, maxiter_Dissipation+1
        !    write (io_lun,'(f10.5)') c(i)
        !  enddo
        !endif
        !db

        allocate (Fdiss(mat_p(matXL(1))%length,nspin), STAT=stat_alloc)
        if (stat_alloc.NE.0) call cq_abort('Error allocating dissipation force: ')
        Fdiss = zero
        if (MDiter.GT.maxiter_Dissipation) call calc_dissipative_force(MDiter,Fdiss)
        !ORI if (MDiter.GE.maxiter_Dissipation) call calc_dissipative_force(MDiter,Fdiss)
        kappa_wk = kappa_diss
        do ispin = 1, nspin
          do len = 1, mat_p(matXL(ispin))%length
            mat_p(matPot(ispin))%matrix(len) = mat_p(matA(ispin))%matrix(len) -  &
                                               mat_p(matXL(ispin))%matrix(len)
            acc = mat_p(matA(ispin))%matrix(len) - mat_p(matXL(ispin))%matrix(len)
            acc_diss = Fdiss(len,ispin)*half/dt
            mat_p(matXLvel(ispin))%matrix(len) = mat_p(matXLvel(ispin))%matrix(len) &
                                                + kappa_wk*half*acc/dt            &
                                                + acc_diss
            mat_p(matXL(ispin))%matrix(len) = mat_p(matXL(ispin))%matrix(len)       + &
                                             dt*mat_p(matXLvel(ispin))%matrix(len) + &
                                             half*kappa_wk*acc                    + &
                                             half*Fdiss(len,ispin)
            mat_p(matXLvel(ispin))%matrix(len) = mat_p(matXLvel(ispin))%matrix(len) &
                                                + kappa_wk*half*acc/dt            &
                                                + acc_diss
          enddo
        enddo
        deallocate (Fdiss, STAT=stat_alloc)
        if (stat_alloc.NE.0) call cq_abort('Error deallocating Fdiss: ')
      else
        kappa_wk = kappa
        do ispin = 1, nspin
          do len = 1, mat_p(matXL(ispin))%length
            mat_p(matPot(ispin))%matrix(len) = mat_p(matA(ispin))%matrix(len) - mat_p(matXL(ispin))%matrix(len)
            acc = mat_p(matA(ispin))%matrix(len) - mat_p(matXL(ispin))%matrix(len)
            mat_p(matXLvel(ispin))%matrix(len) = mat_p(matXLvel(ispin))%matrix(len) &
                                                + kappa_wk*half*acc/dt
            mat_p(matXL(ispin))%matrix(len) = mat_p(matXL(ispin))%matrix(len)       + &
                                             dt*mat_p(matXLvel(ispin))%matrix(len) + &
                                             half*kappa_wk*acc
            mat_p(matXLvel(ispin))%matrix(len) = mat_p(matXLvel(ispin))%matrix(len) &
                                                + kappa_wk*half*acc/dt
          enddo
        enddo
      endif
      ! Calculate kinetic term of electronic degrees of freedom    ! 01/10/2013
      Pot_el =  half*matrix_product_trace(matPot(1) ,matPot(1) )
      KE_el  =  half*matrix_product_trace(matXLvel(1),matXLvel(1))
      if (myid.EQ.0) write (io_lun,'(a,2x,i8,2x,a,f25.15,2x,a,f25.15)') &
                            "    *** MD iter",MDiter,"Tr[(dX/dt)^2]:",KE_el,"Tr[(Z-X)^2]", &
                            Pot_el

      ! Free matPot
      do ispin = nspin, 1, -1
        call free_temp_matrix(matPot(ispin))
      enddo

      if (myid.EQ.0 .AND. iprint_MD.GT.1) &
        write (io_lun,*) "X-matrix propagated via velocity Verlet"

      return
    end subroutine vVerlet_matXL
    !!***

  end module XLBOMD_module
