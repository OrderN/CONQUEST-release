!-*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! -----------------------------------------------------------------------------
! $Id: DiagModule.f90,v 1.9.2.3 2006/03/31 12:04:27 drb Exp $
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
!!   This whole module is discussed in exhaustive detail in the Conquest notes 
!! "Implementation of Diagonalisation in Conquest" (Diagonalisation.tex)
!!  USES
!!   common, cover_module, datatypes, fdf, GenBlas, GenComms, global_module, group_module, matrix_data, matrix_module, 
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
!!    Many changes over the last month - the code now works for 8 and 64
!!    atoms on 1,2,3 and 4 processors in every grid combination I can think
!!    of.  Main changes were to PrepareSend and DistributeCQtoSC.  There will
!!    be more changes as k-point strategy is implemented and as we build the
!!    K matrix.  Also generalised to work row-by-row (i.e. not making blocks 
!!    and atoms commensurate).
!!   16/04/2002 dave
!!    Added distance between atoms to loc, and added k-point reading and loops to FindEvals
!!   16/04/2002 dave
!!    Shifted reading out of FindEvals into readDiagInfo
!!   23/04/2002 dave
!!    Started creation of occupancy and K matrix building routines
!!   23/04/2002 dave
!!    Defined a new derived type to hold arrays relating to distribution - this will allow 
!!    easy incorporation into Conquest
!!   01/05/2002 drb 
!!    Imported to Conquest
!!   17/06/2002 dave
!!    Fixed bugs, tidied (moved initialisation into separate routine etc), added logical 
!!    diagon for solution method choice (maybe should be elsewhere !)
!!   31/07/2002 dave
!!    Added calculation for the matrix M12 which makes the Pulay force (contribution to
!!    dE/dphi_i for the variation of S)
!!   15:37, 03/02/2003 drb 
!!    Various bits of tidying, increase of precision and general improvement
!!   08:20, 2003/07/28 dave
!!    Fairly major reworking to reflect need for multiple CQ elements to fold into one SC element.
!!    New derived type (element), renamed and reworked parts of DistributeData type, reworked PrepareSend
!!    and DistributeCQ_to_SC to reflect all this.
!!   11:48, 30/09/2003 drb 
!!    Changed iprint levels throughout
!!***
module DiagModule

  use datatypes
  use GenComms, ONLY: cq_abort
  use numbers, ONLY: zero

  implicit none

  save 

!!****s* DiagModule/location *
!!  NAME
!!   location
!!  PURPOSE
!!   Holds the location of a matrix element in data_H.  When converting between
!!   Conquest storage and Scalapack storage, we need to know where in data_H a
!!   given row and column of the matrix are, along with support function positions
!!   in that block.  This is stored in loc.
!! 
!!   Now holds the distance between the atoms (true distance, not FSC) for general 
!!   use with k-points (16/04/2002 dave)
!!  AUTHOR
!!   D.R.Bowler
!!  SOURCE
!!
  type location
     integer :: supfn_row, supfn_col
     integer :: loci
     integer :: locj
     real(double) :: dx,dy,dz
  end type location
!!***

!!****s* DiagModule/elements *
!!  NAME
!!   elements
!!  PURPOSE
!!   Holds an array of locations of matrix elements.  When redistributing Conquest matrix elements
!!   to Scalapack matrix elements, this will be required when the support function radius is more 
!!   than quarter of the unit cell side (i.e. when an atom and its periodic image are both neighbours
!!   of the same primary set atom).  Allows more than one CQ element to contribute to a given SC 
!!   element
!!  AUTHOR
!!   D.R.Bowler
!!  SOURCE
!!
  type element
     integer :: n_elements
     type(location), dimension(:), pointer :: where
  end type element
!!***

!!****s* DiagModule/DistributeData *
!!  NAME
!!   DistributeData
!!  PURPOSE
!!   Holds arrays related to distribution of a matrix from Conquest compressed row storage to
!!   ScaLAPACK 2D block cyclic storage
!!  AUTHOR
!!   D.R.Bowler
!!  SOURCE
!!
  type DistributeData
     integer, dimension(:), pointer :: num_rows, start_row ! Receiving data
     type(element), dimension(:,:,:), pointer :: images
     integer, dimension(:), pointer :: send_rows, firstrow
  end type DistributeData
!!***

!!****s* DiagModule/Krecv_data *
!!  NAME
!!   Krecv_data
!!  PURPOSE
!!   Holds various data about matrix elements we're going to build from eigenvectors we'll receive
!!  AUTHOR
!!   D.R.Bowler
!!  SOURCE
!!
    type Krecv_data
       integer, dimension(:), pointer :: ints  ! How many interactions ?
       !integer, dimension(:), pointer :: ndimi  ! Dimension of i
       integer, dimension(:), pointer :: ndimj  ! Dimension of j
       integer :: orbs  ! How many orbitals ?
       integer, dimension(:,:), pointer :: prim_atom ! Primary set atom for a given interaction
       integer, dimension(:,:), pointer :: locj ! Where to put matrix element in data_Matrix
       real(double), dimension(:,:), pointer :: dx, dy, dz ! Vector between atoms
    end type Krecv_data
!!***

  ! The matrix that holds the SC data - the HAMILTONIAN 
  complex(double_cplx), dimension(:,:), allocatable :: SCHmat, SCSmat
  complex(double_cplx), dimension(:,:), allocatable :: z
  ! Buffer for receiving data
  complex(double_cplx), dimension(:,:), allocatable :: RecvBuffer
  ! Buffer for sending data
  complex(double_cplx), dimension(:,:), allocatable :: SendBuffer
  ! Data type that stores details of where to put and where to get elements
  type(DistributeData) :: DistribH, DistribS

  ! K-point data - here so that reading of k-points can take place in different routine to FindEvals
  integer :: nkp
  real(double), allocatable, dimension(:,:) :: kk
  real(double), allocatable, dimension(:) :: wtk
  real(double), allocatable, dimension(:,:) :: occ
  real(double), parameter :: kT = 0.0001_double

  logical :: first = .true.

  logical :: diagon ! Do we diagonalise or use O(N) ?
  ! -------------------------------------------------------
  ! RCS ident string for object file id
  ! -------------------------------------------------------
  character(len=80), private :: RCSid = "$Id: DiagModule.f90,v 1.9.2.3 2006/03/31 12:04:27 drb Exp $"

  ! Local scratch data
  real(double), allocatable, dimension(:,:) :: w
  !complex(double_cplx), dimension(:),allocatable :: work, rwork, gap
  complex(double_cplx), dimension(:),allocatable :: work
  real(double), dimension(:),allocatable :: rwork, gap
  integer, dimension(:), allocatable :: iwork, ifail, iclustr

  ! BLACS variables
  integer :: me, context

contains

! -----------------------------------------------------------------------------
! Subroutine FindEvals
! -----------------------------------------------------------------------------

!!****f* DiagModule/FindEvals *
!!
!!  NAME 
!!   FindEvals - finds the eigenvalues
!!  USAGE
!! 
!!  PURPOSE
!!   Call the ScalapackFormat routines and DiagModule routines to 
!!   find the eigenvalues by exact diagonalisation.  See the Conquest notes 
!! "Implementation of Diagonalisation within Conquest" for a detailed 
!!   discussion of this routine.
!!
!!   Note added 31/07/2002 on Pulay force (drb)
!! 
!!   The Pulay contribution to the force (which is also the rate of change of the total energy with respect
!!   to a given support function) requires the matrix M12 (which has range S).  When using LNV order N this
!!   is given by 3LHL - 2(LHLSL + LSLHL), but during diagonalisation it's given by:
!!
!!   M12_ij = -\sum_k w_k \sum_n f_n \epsilon_n c^n_i c^n_j
!!
!!   This can be trivially implemented using buildK if we scale the occupancies by the eigenvalues after 
!!   building K and before building M12, and this is what is done.  Within buildK, the c_j coefficients are
!!   scaled by f_n before a dot product is taken with c_i, so scaling the occupancies by the eigenvalues 
!!   effectively builds M12.
!!  INPUTS
!!   real(double) :: electrons - number of electrons in system
!!  USES
!!   common, datatypes, GenComms, global_module, matrix_data, maxima_module, primary_module, ScalapackFormat
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
!!    Added reading of k-points (using fdf) and loop over k-points for the
!!    diagonalisation calls.  Also now passes k-point to DistributeCQtoSC and calls a new
!!    ScaLAPACK routine (pzheev) to do the diagonalisation (as the Hamiltonian is now complex)
!!   16/04/2002 dave
!!    Moved all I/O and allocation for k-points down to readDiagInfo
!!   18/04/2002 dave
!!    Changed DistributeCQ_to_SC so that the matrix we're distributing into is passed - this will
!!    allow distribution of S as well as H
!!   22/04/2002 dave
!!    Tidied and moved code around in preparation for building K from eigenvectors
!!   23/04/2002 dave
!!    Added calls to find Fermi level and occupancies of k points, and moved output of eigenvalues
!!    so that occupancies can be written out as well
!!   01/05/2002 dave
!!    Moved allocation of SCHmat and z out of PrepareRecv
!!   29/05/2002 dave
!!    Moved initialisation (BLACS start up, descinit, DistribH, allocation of memory) to initDiag
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
!!    Removed all explicit references to data_ variables and rewrote in terms of new 
!!    matrix routines
!!  SOURCE
!!
  subroutine FindEvals(electrons)

    use datatypes
    use numbers, ONLY: zero
    use units
    use global_module, ONLY: iprint_DM, ni_in_cell, numprocs, area_DM
    use GenComms, ONLY: my_barrier, cq_abort, mtmini, mtime, gcopy, myid
    use ScalapackFormat, ONLY: matrix_size, proc_rows, proc_cols, deallocate_arrays, &
         block_size_r, block_size_c, my_row
    use mult_module, ONLY: matH, matS, matK, matM12, matrix_scale, matrix_product_trace
    use matrix_data, ONLY: Hrange, Srange
    use primary_module, ONLY: bundle
    use species_module, ONLY: species, nsf_species
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl, type_int

    implicit none

    ! Passed variables
    real(double) :: electrons
    
    ! Local variables
    real(double) :: Ef, bandE, abstol, a
    real(double), external :: dlamch
    complex(double_cplx), dimension(:,:), allocatable :: expH
    integer :: merow, mecol, info, lwork, stat, row_size, nump
    integer :: iunit, i, j, k, l, lrwork, nsf1, k1, k2, col_size, np
    integer :: liwork, m, mz, prim_size
    integer, dimension(50) :: desca,descz,descb

    if(iprint_DM>=2.AND.myid==0) write(*,fmt='(10x,"Entering FindEvals")')
    ! Read appropriate data for Scalapack diagonalisation - k-points, block sizes etc
    !if(first) call readDiagInfo(proc_rows,proc_cols,block_size_r,block_size_c)
    matrix_size = 0 !ni_in_cell*nsf
    do i=1,ni_in_cell
       matrix_size = matrix_size + nsf_species(species(i))
    end do
    prim_size = 0
    do i=1,bundle%n_prim
       prim_size = prim_size + nsf_species(bundle%species(i))
    end do
    ! Check for block size factoring into matrix size
    a = real(matrix_size)/real(block_size_r)
    if(a - real(floor(a))>1e-8_double) call cq_abort('block_size_r not a factor of matrix size ! ',matrix_size, block_size_r)
    a = real(matrix_size)/real(block_size_c)
    if(a - real(floor(a))>1e-8_double) call cq_abort('block_size_c not a factor of matrix size ! ',matrix_size, block_size_c)
    ! Initialise - start BLACS, sort out matrices, allocate memory
    call initDiag(context, me, desca, descb, descz, lwork, lrwork, liwork)
    !my_row = 0
    ! -------------------------------------------------------------------------------------------------
    ! Start diagonalisation
    ! -------------------------------------------------------------------------------------------------
    ! First diagonalisation - get eigenvalues only (so that we can find Ef)
    call mtmini
    abstol = 1e-30_double!pdlamch(context,'U')
    if(iprint_DM>=2.AND.myid==0) write(*,fmt='(10x,"In FindEvals, tolerance is ",g20.12)') abstol
    do i=1,nkp ! Build and send Hamiltonian for each k-point
       ! Form the Hamiltonian for this k-point and send it to appropriate processors
       if(iprint_DM>=3) write(*,*) myid,' Calling DistributeCQ_to_SC for H'
       call my_barrier()
       call DistributeCQ_to_SC(DistribH,matH,kk(1,i),SCHmat)
       ! Form the overlap for this k-point and send it to appropriate processors
       if(iprint_DM>=3) write(*,*) myid,' Calling DistributeCQ_to_SC for S'
       call DistributeCQ_to_SC(DistribS,matS,kk(1,i),SCSmat)
       if(iprint_DM>=3) write(*,*) myid,' Done DistributeCQ_to_SC for S'
       ! Now, if this processor is involved, do the diagonalisation
       if(iprint_DM>=3) write(*,*) myid,'Proc row, cols, me: ',proc_rows, proc_cols, me, i, nkp
       if(me<proc_rows*proc_cols) then
          ! Call the diagonalisation routine for generalised problem H.psi = E.S.psi
          call pzhegvx(1,'N','A','U',matrix_size,SCHmat,1,1,desca,SCSmat,1,1,descb,&
               0.0d0,0.0d0,0,0,abstol,m,mz,w(1,i),-1.0_double,z,1,1,descz,&
               work,lwork,rwork,lrwork,iwork,liwork,ifail,iclustr,gap,info)
          !call zhegvx(1,'N','A','U',matrix_size,SCHmat,matrix_size,SCSmat,matrix_size,&               
          !     0.0d0,0.0d0,0,0,abstol,m,w(1,i),z,matrix_size,&
          !     work,lwork,rwork,iwork,ifail,info)
          if(info/=0) call cq_abort("FindEvals: pzheev failed !",info)
       end if ! if(me<proc_rows*proc_cols)
    end do ! End do i=1,nkp
    if(iprint_DM>=2) write(*,2) myid,mtime()
    ! Find Fermi level, given the eigenvalues at all k-points (in w)
    if(me<proc_rows*proc_cols) then
       call findFermi(electrons,w,matrix_size,nkp,Ef)
       call gcopy(Ef)
       call gcopy(occ,matrix_size,nkp)
    else
       call gcopy(Ef)
       call gcopy(occ,matrix_size,nkp)
    end if
    ! Now write out eigenvalues and occupancies
    if(iprint_DM>=3.AND.myid==0) then
       bandE = 0.0_double
       do i=1,nkp
          write(*,7) i,kk(1,i),kk(2,i),kk(3,i)
          do j=1,matrix_size,3
             if (j==matrix_size) then
                write(*,8) w(j,i),occ(j,i)
                bandE = bandE + w(j,i)*occ(j,i)
             else if (j==matrix_size-1) then
                write(*,9) w(j,i),occ(j,i),w(j+1,i),occ(j+1,i)
                bandE = bandE + w(j,i)*occ(j,i)+ w(j+1,i)*occ(j+1,i)
             else
                write(*,10) w(j,i),occ(j,i),w(j+1,i),occ(j+1,i),w(j+2,i),occ(j+2,i)
                bandE = bandE + w(j,i)*occ(j,i)+ w(j+1,i)*occ(j+1,i)+ w(j+2,i)*occ(j+2,i)
             endif
          end do ! j=matrix_size
       end do ! End do i=1,nkp
       write(*,4) bandE
    end if ! if(iprint_DM>=1.AND.myid==0)
    ! Allocate space to expand eigenvectors into (i.e. when reversing ScaLAPACK distribution)
    allocate(expH(matrix_size,prim_size),STAT=stat)
    if(stat/=0) call cq_abort('FindEvals: failed to alloc expH',matrix_size)
    call reg_alloc_mem(area_DM,2*matrix_size*prim_size,type_dbl)
    call mtmini
    call matrix_scale(zero,matK)
    call matrix_scale(zero,matM12)
    ! Second diagonalisation - get eigenvectors and build K
    do i=1,nkp
       ! Form the Hamiltonian for this k-point and send it to appropriate processors
       if(iprint_DM>=3) write(*,*) myid,' Calling DistributeCQ_to_SC for H'
       call DistributeCQ_to_SC(DistribH,matH,kk(1,i),SCHmat)
       ! Form the overlap for this k-point and send it to appropriate processors
       if(iprint_DM>=3) write(*,*) myid,' Calling DistributeCQ_to_SC for S'
       call DistributeCQ_to_SC(DistribS,matS,kk(1,i),SCSmat)
       ! Now, if this processor is involved, do the diagonalisation
       if(me<proc_rows*proc_cols) then
          ! Call the diagonalisation routine for generalised problem H.psi = E.S.psi
          call pzhegvx(1,'V','A','U',matrix_size,SCHmat,1,1,desca,SCSmat,1,1,descb,&
               0.0d0,0.0d0,0,0,abstol,m,mz,W(1,i),-1.0_double,z,1,1,descz,&
               work,lwork,rwork,lrwork,iwork,liwork,ifail,iclustr,gap,info)
          !call zhegvx(1,'V','A','U',matrix_size,SCHmat,matrix_size,SCSmat,matrix_size,&
          !     0.0d0,0.0d0,0,0,abstol,m,w(1,i),z,matrix_size,&
          !     work,lwork,rwork,iwork,ifail,info)
          if(info<0) call cq_abort("FindEvals: pzheev failed !",info)
          if(info>=1) write(*,*) 'Problem - info returned as: ',info
       end if ! if(me<proc_rows*proc_cols)
       if(iprint_DM>=5) write(*,*) myid,' Calling barrier'
       call my_barrier()
       if(iprint_DM>=5) write(*,*) myid,' Calling DistributeSC_to_Ref'
       ! Reverse the CQ to SC distribution so that eigenvector coefficients for atoms are on the appropriate processor
       call DistributeSC_to_ref(DistribH,z,expH) 
       ! Build K from the eigenvectors
       if(iprint_DM>=4) write(*,*) myid,' Calling buildK ',Hrange, matK
       call buildK(Hrange, matK, occ(1:matrix_size,i), kk(1:3,i), wtk(i), expH)
       ! Build matrix needed for Pulay force
       ! We scale the occupation number for this k-point by the eigenvalues in order to build the matrix M12
       ! We can do this simply because we won't use them again (though we could use a dummy variable if we
       ! wanted to use them again)
       do j=1,matrix_size 
          occ(j,i) = -occ(j,i)*w(j,i)
       end do
       ! Now build data_M12_ij (=-\sum_n eps^n c^n_i c^n_j - hence scaling occs by eps allows reuse of buildK)
       call buildK(Srange, matM12, occ(1:matrix_size,i), kk(1:3,i), wtk(i), expH)
    end do ! End do i=1,nkp
    if(iprint_DM>=2) write(*,3) myid,mtime()
    ! -------------------------------------------------------------------------------------------------
    ! End diagonalisation
    ! -------------------------------------------------------------------------------------------------
    ! Write out the band energy and trace of K
    if(iprint_DM>=1) then
       bandE = 2.0_double*matrix_product_trace(matK,matH)
       if(myid==0) write(*,5) en_conv*bandE,en_units(energy_units)
       bandE = 2.0_double*matrix_product_trace(matS,matM12)
       if(myid==0) write(*,6) en_conv*bandE,en_units(energy_units)
    end if
    call my_barrier()
    ! Deallocate space
    deallocate(expH,z,SCSmat,SCHmat,STAT=stat)
    if(stat/=0) call cq_abort('Error deallocating (0a) in FindEvals',stat)
    call reg_dealloc_mem(area_DM,2*matrix_size*prim_size+6*row_size*col_size,type_dbl)
    deallocate(w,occ,STAT=stat)
    if(stat/=0) call cq_abort('Error deallocating (0b) in FindEvals',stat)
    call reg_dealloc_mem(area_DM,2*matrix_size*nkp,type_dbl)
    ! Shut down BLACS
    if(me<proc_rows*proc_cols) then
       call blacs_gridexit(context)
       deallocate(work,rwork,iwork,ifail,iclustr,gap,STAT=stat)
       call reg_dealloc_mem(area_DM,2*lwork+lrwork+liwork,type_dbl)
       call reg_dealloc_mem(area_DM,matrix_size+2*proc_rows*proc_cols+proc_rows*proc_cols,type_int)
       !deallocate(work,rwork,iwork,ifail,STAT=stat)
       if(stat/=0) call cq_abort('Error deallocating (1) in FindEvals',stat)
    end if
    do i=1,size(DistribS%images,3)
       do j=1,size(DistribS%images,2)
          do k=1,size(DistribS%images,1)
             if(DistribS%images(k,j,i)%n_elements>0) then
                deallocate(DistribS%images(k,j,i)%where,STAT=stat)
                if(stat/=0) call cq_abort('Error deallocating (2) in FindEvals',stat)
             end if
          end do
       end do
    end do
    deallocate(DistribS%images,STAT=stat)
    if(stat/=0) call cq_abort('Error deallocating (2) in FindEvals',stat)
    deallocate(DistribS%num_rows,DistribS%start_row,DistribS%send_rows,DistribS%firstrow,STAT=stat)
    if(stat/=0) call cq_abort('Error deallocating (2a) in FindEvals',stat)
    call reg_dealloc_mem(area_DM,4*numprocs,type_int)
    do i=1,size(DistribH%images,3)
       do j=1,size(DistribH%images,2)
          do k=1,size(DistribH%images,1)
             if(DistribH%images(k,j,i)%n_elements>0) then
                deallocate(DistribH%images(k,j,i)%where,STAT=stat)
                if(stat/=0) call cq_abort('Error deallocating (2) in FindEvals',stat)
             end if
          end do
       end do
    end do
    deallocate(DistribH%images,STAT=stat)
    if(stat/=0) call cq_abort('Error deallocating (2) in FindEvals',stat)
    deallocate(DistribH%num_rows,DistribH%start_row,DistribH%send_rows,DistribH%firstrow,STAT=stat)
    if(stat/=0) call cq_abort('Error deallocating (2a) in FindEvals',stat)
    call deallocate_arrays
    return
2   format(10x,'Proc: ',i5,' Time taken for eval diag: ',f15.8,' ms')
3   format(10x,'Proc: ',i5,' Time taken for evec diag: ',f15.8,' ms')
4   format(10x,'Sum of eigenvalues: ',f18.11,' ',a2)
5   format(10x,'Energy as 2Tr[K.H]: ',f18.11,' ',a2)
6   format(10x,'2Tr[S.G]: ',f18.11,' ',a2)
7   format(10x,'Eigenvalues and occupancies for k-point ',i3,' : ',3f12.5)
8   format(10x,f12.5,f6.3,2x)
9   format(10x,f12.5,f6.3,2x,f12.5,f6.3,2x)
10  format(10x,f12.5,f6.3,2x,f12.5,f6.3,2x,f12.5,f6.3,2x)
11  format(10x,'Proc: ',i5,' Diagonalising for eigenvectors')
12  format(10x,'Proc: ',i5,' row, col size: ',2i5)
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
!!   call initDiag(context, me, desca, descb, descz, work, rwork, gap, iwork, ifail, iclustr, lwork, lrwork, liwork)
!!  PURPOSE
!!   Contains various routines and assignments that initialise the diagonalisation
!!  INPUTS
!!   integer :: context - BLACS context handle
!!   integer :: me      - BLACS processor id
!!   integer, dimension(50) :: desca, descb, descz - descriptors for H, S and eigenvectors
!!
!!   All variables below are ScaLAPACK scratch space and scratch sizes (found and allocated here)
!!
!!   complex(double_cplx), dimension(:),allocatable :: work, rwork, gap
!!   integer, dimension(:), allocatable :: iwork, ifail, iclustr
!!   integer :: lwork, lrwork, liwork
!!  USES
!!   ScalapackFormat, matrix_data, GenComms
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   29/05/2002
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine initDiag(context, me, desca, descb, descz, lwork, lrwork, liwork)

    use ScalapackFormat, ONLY: allocate_arrays,ref_to_SC_blocks,make_maps,find_SC_row_atoms, &
         find_ref_row_atoms,find_SC_col_atoms, proc_start, block_size_r, block_size_c, proc_rows, proc_cols, &
         matrix_size
    use global_module, ONLY: iprint_DM, numprocs
    use matrix_data, ONLY: Hrange, Srange, mat
    use GenComms, ONLY: my_barrier, cq_abort, myid
    use mult_module, ONLY: matH, matS

    implicit none

    ! Passed variables
    integer :: context, me
    integer, dimension(50) :: desca, descb, descz
    integer :: lwork, lrwork, liwork

    ! Local variables
    integer :: row_size, col_size, nump, merow, mecol, info, stat, m, mz
    real(double) :: rwo(1)
    complex(double) :: wo(1)
    integer :: iwo(1)

    ! ScalapackFormat - do various translation tasks between Conquest compressed row
    ! storage and distributed Scalapack form. The arrays created will be used to form
    ! the Hamiltonian and K matrices.
    call allocate_arrays
    call ref_to_SC_blocks
    call make_maps
    call find_SC_row_atoms
    call find_ref_row_atoms
    call find_SC_col_atoms
    ! Now that's done, we can prepare to distribute things
    call my_barrier()
    call PrepareRecv(DistribH)
    call PrepareSend(matH,Hrange,DistribH)
    call my_barrier()
    call PrepareRecv(DistribS)
    call PrepareSend(matS,Srange,DistribS)
    call my_barrier()
    ! First, work out how much data we're going to receive
    ! How many rows and columns do we have ?
    row_size = proc_start(myid+1)%rows*block_size_r ! Sizes of local "chunk"
    col_size = proc_start(myid+1)%cols*block_size_c
    if(iprint_DM>=3) write(*,12) myid,row_size,col_size
    ! Allocate space for the distributed Scalapack matrices
    allocate(SCHmat(row_size,col_size),STAT=stat)
    if(stat/=0) call cq_abort("DiagModule: Could not alloc SCHmat",stat)
    allocate(SCSmat(row_size,col_size),STAT=stat)
    if(stat/=0) call cq_abort("DiagModule: Could not alloc SCSmat",stat)
    allocate(z(row_size,col_size),STAT=stat)
    if(stat/=0) call cq_abort("DiagModule: Could not alloc z",stat)
    SCHmat = 0.0_double
    SCSmat = 0.0_double
    z = 0.0_double
    ! Allocate eigenvalue storage
    allocate(w(matrix_size,nkp),occ(matrix_size,nkp),STAT=stat)
    if(stat/=0) call cq_abort('FindEvals: failed to alloc w and occ',matrix_size,nkp)
    ! Start up BLACS
    nump = numprocs
    me = myid
    call blacs_pinfo(me,nump)
    call blacs_get(-1,0,context)
    call blacs_gridinit(context,'R',proc_rows,proc_cols)
    call blacs_gridinfo(context,proc_rows,proc_cols,merow,mecol)
    if(iprint_DM>=3) write(*,1) myid,me,merow,mecol
    call blacs_pinfo(me,nump)
    call my_barrier
    if(nump/=numprocs) call cq_abort('FindEvals: BLACS procs error !',nump,numprocs)
    ! If participating in diagonalisation, register matrices and find work space requirements
    if(me<proc_rows*proc_cols) then
       ! Register the description of the distribution of H
       call descinit(desca,matrix_size,matrix_size,&
            block_size_r,block_size_c,&
            0,0,context,row_size,info)
       if(info/=0) call cq_abort("FindEvals: descinit(a) failed !",info)
       ! Register the description of the distribution of S
       call descinit(descb,matrix_size,matrix_size,&
            block_size_r,block_size_c,&
            0,0,context,row_size,info)
       if(info/=0) call cq_abort("FindEvals: descinit(a) failed !",info)
       ! And register eigenvector distribution
       call descinit(descz,matrix_size,matrix_size,&
            block_size_r,block_size_c,&
            0,0,context,row_size,info)
       ! Find scratch space requirements for ScaLAPACk
       if(info/=0) call cq_abort("FindEvals: descinit(z) failed !",info)
       allocate(ifail(matrix_size),iclustr(2*proc_rows*proc_cols), gap(proc_rows*proc_cols),STAT=stat)
       !allocate(work(1),rwork(1),iwork(1),STAT=stat)
       if(stat/=0) call cq_abort('FindEvals: failed to alloc work and rwork',lwork,lrwork)
       call pzhegvx(1,'V','A','U',matrix_size,SCHmat,1,1,desca,SCSmat,1,1,descb,&
            0.0d0,0.0d0,0,0,1.0e-307_double,m,mz,w(1,1),-1.0_double,z,1,1,descz,&
            wo,-1,rwo,-1,iwo,-1,ifail,iclustr,gap,info)
       ! Allocate scratch space for ScaLAPACK
       lwork  = 2*real(wo(1))
       lrwork = 2*rwo(1)
       liwork = iwo(1)
       !deallocate(work,rwork,iwork,STAT=stat)
       !if(stat/=0) call cq_abort('FindEvals: failed to alloc work and rwork',lwork,lrwork)
       !lwork = 9*matrix_size
       stat = 0
       allocate(work(lwork),STAT=stat)
       if(stat/=0) call cq_abort('FindEvals: failed to alloc work',lwork)
       allocate(rwork(lrwork),STAT=stat)
       if(stat/=0) call cq_abort('FindEvals: failed to alloc rwork',lrwork)
       allocate(iwork(liwork),STAT=stat)
       !allocate(work(lwork),rwork(7*matrix_size),iwork(5*matrix_size),ifail(matrix_size),STAT=stat)
       if(stat/=0) call cq_abort('FindEvals: failed to alloc iwork',liwork)
    end if
    call my_barrier
    return
1   format(10x,'Proc: ',i5,' BLACS proc, row, col: ',3i5)
12  format(10x,'Proc: ',i5,' row, col size: ',2i5)
  end subroutine initDiag
!!***

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
!!    Changed to use the DistributeData derived type which is passed down. Also removed max_rows.
!!   2004/11/10 drb
!!    Changed nsf to come from maxima, not common
!!   2006/08/30 16:49 dave
!!    Added allocate for arrays in Distrib
!!  SOURCE
!!
  subroutine PrepareRecv(Distrib)

    use ScalapackFormat, ONLY: proc_start, SC_row_block_atom, &
         block_size_r, block_size_c, blocks_r, my_row, matrix_size
    use global_module, ONLY: iprint_DM, ni_in_cell, numprocs
    use group_module, ONLY: parts
    use GenComms, ONLY: myid, cq_abort

    implicit none

    ! Passed variables
    type(DistributeData), intent(out) :: Distrib

    ! Local variables
    integer :: i,j,stat
    integer :: row_size,col_size,count
    integer :: row, rowblock,proc

    if(iprint_DM>=2.AND.myid==0) write(*,fmt='(10x,"Entering PrepareRecv")')
    allocate(Distrib%num_rows(numprocs),Distrib%start_row(numprocs),Distrib%send_rows(numprocs),&
         Distrib%firstrow(numprocs),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating Distrib arrays: ",numprocs,stat)
    ! Prepare counters
    Distrib%num_rows = 0
    count = 1
    i = 1
    do rowblock=1,blocks_r
       if(my_row(rowblock)>0) then ! If this row block is part of my chunk
          do row = 1,block_size_r
             if(iprint_DM>=5) write(*,4) myid,i,rowblock,row
             ! Find processor and increment processors and rows from proc
             proc = parts%i_cc2node(SC_row_block_atom(row,rowblock)%part)
             if(Distrib%num_rows(proc)==0) Distrib%start_row(proc)=count ! Where data from proc goes
             Distrib%num_rows(proc) = Distrib%num_rows(proc)+1
             if(iprint_DM>=5) write(*,3) myid,i,proc,Distrib%start_row(proc),Distrib%num_rows(proc)
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
!!   Sets up arrays which are vital to the sending of data when going from
!!   CQ format (compressed arrays) to SC format.  Works in three stages:
!!
!!   i) Loops over SC matrix rows by going over atoms and support functions.  Notes
!!      which processors each row is sent to, and where the data for a processor starts.
!!   ii) Works out how many rows and columns are sent, and allocates data storage
!!   iii) For each non-zero matrix element, works out where in SC format matrix it goes
!!  INPUTS
!!   type(matrix) :: mat
!!   type(DistributeData) :: Distrib
!!  USES
!!   global_module, common, maxima_module, matrix_module, group_module, primary_module,
!!   cover_module, ScalapackFormat
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
!!    Updated the early stages where we're trying to work out which processors
!!    we send each row to, and where the data for each processor starts.  Simplified
!!    the algorithm significantly. Also added further comments and changed matrix from
!!    Hmat to mat.
!!   16/04/2002 dave
!!    Added distance vector between atoms to loc so that the phases to be applied to the
!!    Hamiltonian can be calculated later
!!   16/04/2002 dave
!!    Added welcome statement
!!   23/04/2002 dave
!!    Added intent to mat
!!   23/04/2002 dave
!!    Added a derived datatype (Distrib) which is passed down to store all the arrays
!!    related to communications - it will allow easier generalisation to Conquest
!!   01/05/2002 dave
!!    Added level 4 iprint_DM to a write statement
!!   08:32, 2003/07/28 dave
!!    Major reworking to allow multiple CQ elements to contribute to a single SC elements (should 
!!    have been there from the start but overlooked !)
!!   2004/10/29 drb
!!    Added check on size of n_elements before allocating Distrib
!!   2004/11/10 drb
!!    Changed nsf to come from maxima, not common
!!   09:09, 11/05/2005 dave 
!!    Added check on j,k dimensions
!!  SOURCE
!!
  subroutine PrepareSend(matA,range,Distrib)

    use global_module, ONLY: numprocs, iprint_DM, x_atom_cell, y_atom_cell, z_atom_cell
    use GenComms, ONLY: myid, my_barrier
    use maxima_module, ONLY: maxnsf
    use matrix_module, ONLY: matrix, matrix_halo
    use group_module, ONLY: parts
    use primary_module, ONLY: bundle
    use cover_module, ONLY: BCS_parts
    use ScalapackFormat, ONLY : CC_to_SC,maxrow,maxcol,proc_block, &
         SC_to_refx,SC_to_refy, block_size_r, block_size_c, blocks_c, proc_start
    use matrix_data, ONLY: mat, halo
    use species_module, ONLY: nsf_species

    implicit none

    ! Passed variables
    integer :: range, matA
    type(DistributeData), intent(out) :: Distrib

    ! Local variables
    integer, allocatable, dimension(:,:,:) :: sendlist
    integer, allocatable, dimension(:,:,:) :: ele_list
    integer, allocatable, dimension(:) :: firstcol, sendto
    integer :: part, memb, neigh, ist, atom_num
    integer :: Row_FSC_part, Row_FSC_seq, Col_FSC_part, Col_FSC_seq
    integer :: SCblockr, SCblockc, SCrowc, i, j, k, l, proc, stat, supfn_r, supfn_c
    integer :: maxr, maxc, CC, row, currow, start, gcspart
    integer :: FSCpart, i_acc_prim_nsf, prim_nsf_so_far

    if(iprint_DM>=2.AND.myid==0) write(*,fmt='(10x,"Entering PrepareSend")')
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
    if(iprint_DM>=2.AND.myid==0) write(*,*) 'Part i ',myid
    call my_barrier()
    ! Work out which processors I send to
    i_acc_prim_nsf = 0
    do part = 1,bundle%groups_on_node ! Loop over primary set partitions
       if(bundle%nm_nodgroup(part)>0) then
          CC = parts%ngnode(parts%inode_beg(myid+1)+part-1)
          do memb = 1,bundle%nm_nodgroup(part) ! Loop over atoms
             atom_num = bundle%nm_nodbeg(part)+memb-1
             Row_FSC_part = CC
             Row_FSC_seq  = memb
             do supfn_r = 1,nsf_species(bundle%species(atom_num))
                i_acc_prim_nsf = i_acc_prim_nsf + 1
                ! Here we work out which row of the unpacked SC matrix we're in
                SCblockr = CC_to_SC(Row_FSC_part,Row_FSC_seq,supfn_r)%block_r
                ! Loop over columns of blocks in SC format
                do i=1,blocks_c
                   proc = proc_block(SC_to_refx(SCblockr,i), SC_to_refy(SCblockr,i))
                   sendlist(supfn_r,atom_num,proc) = 1
                   ! If this is the first row and column sent to this processor, then record where data starts
                   if(sendto(proc)==0) then
                      sendto(proc) = 1
                      Distrib%firstrow(proc) = i_acc_prim_nsf
                      firstcol(proc)=(i-1)*block_size_c+1
                   end if ! End if sendto == 0
                end do ! End do i=1,blocks_c
             end do ! End do supfn_r
          end do ! End do memb = 1,bundle%nm_nodgroup
       end if ! End if nm_nodgroup>0
    end do ! End do part = 1,groups_on_node
    ! Part (ii)
    ! Work out how many rows to send to each remote processor
    if(iprint_DM>=2.AND.myid==0) write(*,*) 'Part ii ',myid
    call my_barrier()
    Distrib%send_rows = 0
    do proc=1,numprocs ! Loop over processors
       start = 0
       if(sendto(proc)/=0) then ! If we send to this processor
          do i=1,bundle%n_prim
             do supfn_r = 1,nsf_species(bundle%species(i))
                if(sendlist(supfn_r,i,proc)==1.AND.start==0) then
                   start=1
                end if
                if(sendlist(supfn_r,i,proc)==1) then
                   Distrib%send_rows(proc)=Distrib%send_rows(proc)+1
                   if(iprint_DM>=5) write(*,11) myid,proc,i,supfn_r
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
       write(*,*) myid,' Allocating ele_list: ',numprocs,maxr,maxc
       call my_barrier()
    end if
    stat = 0
    allocate(ele_list(numprocs,maxr,maxc),STAT=stat)
    if(stat/=0) call cq_abort("DiagModule: Failed to alloc ele_list ",stat)
    ele_list = 0
    if(iprint_DM>=2.AND.myid==0) write(*,*) 'Part iia ',myid
    call my_barrier()
    ! Added to count elements: how many CQ elements fall onto this one SC one ?
    !write(*,*) myid,' Parts: ',bundle%groups_on_node
    i_acc_prim_nsf = 0
    do part = 1,bundle%groups_on_node ! Loop over primary set partitions
       if(bundle%nm_nodgroup(part)>0) then
          CC = parts%ngnode(parts%inode_beg(myid+1)+part-1)
          !write(*,*) myid,' Membs: ',part,bundle%nm_nodgroup(part)
          do memb = 1,bundle%nm_nodgroup(part) ! Loop over atoms
             atom_num = bundle%nm_nodbeg(part)+memb-1
             Row_FSC_part = CC
             Row_FSC_seq  = memb
             !write(*,*) myid,' Nabs:  ',part,memb,mat(part,range)%n_nab(memb)
             prim_nsf_so_far = i_acc_prim_nsf
             do neigh = 1, mat(part,range)%n_nab(memb) ! Loop over neighbours of atom
                if(iprint_DM>=5) write(*,3) myid,neigh
                ist = mat(part,range)%i_acc(memb)+neigh-1
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
                      ! Find processor that this block belongs to
                      proc = proc_block(SC_to_refx(SCblockr,SCblockc), &
                           SC_to_refy(SCblockr,SCblockc))
                      ! Create atom numbers and store posn in array at point
                      j = row-Distrib%firstrow(proc)+1  
                      k = (SCblockc-1)*block_size_c + SCrowc-firstcol(proc)+1
                      if(j>maxr.OR.k>maxc.OR.proc>numprocs) then
                         write(*,*) 'Error ! Maxr/c exceeded: ',j,maxr,k,maxc,proc,numprocs
                         call cq_abort('Error in counting elements !')
                      end if
                      if(iprint_DM>=5) then
                         write(*,*) myid,' Loop: ',part,memb,supfn_r,supfn_c
                         write(*,*) myid,' j,k,proc: ',j,k,proc
                      end if
                      ele_list(proc,j,k) = ele_list(proc,j,k) + 1
                   end do
                end do ! End do supfn_r
             end do ! End do neigh
          end do ! End do memb = 1,bundle%nm_nodgroup
       end if ! End if nm_nodgroup>0
    end do ! End do part = 1,groups_on_node
    call my_barrier()
    if(iprint_DM>=5) then
       write(*,10) myid,maxr,maxc,maxrow,maxcol
       call my_barrier()
    end if
    if(iprint_DM>=2.AND.myid==0) write(*,*) 'Part iib ',myid
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
       write(*,*) myid," calling part iii"
       call my_barrier()
    end if
    ! Part (iii)
    ! Build the loc index - for each non-zero element in data_H, record where it goes in SC matrix
    if(iprint_DM>=2.AND.myid==0) write(*,*) 'Part iii ',myid
    call my_barrier()
    i_acc_prim_nsf = 0
    do part = 1,bundle%groups_on_node ! Loop over primary set partitions
       if(iprint_DM>=5) write(*,1) myid,part
       if(bundle%nm_nodgroup(part)>0) then ! If there are atoms in partition
          CC = parts%ngnode(parts%inode_beg(myid+1)+part-1)
          do memb = 1,bundle%nm_nodgroup(part) ! Loop over atoms
             if(iprint_DM>=5) write(*,2) myid,memb
             atom_num = bundle%nm_nodbeg(part)+memb-1
             Row_FSC_part = CC
             Row_FSC_seq  = memb
             prim_nsf_so_far = i_acc_prim_nsf
             do neigh = 1, mat(part,range)%n_nab(memb) ! Loop over neighbours of atom
                if(iprint_DM>=5) write(*,3) myid,neigh
                ist = mat(part,range)%i_acc(memb)+neigh-1
                ! Establish FSC number of neighbour
                Col_FSC_part = BCS_parts%lab_cell(mat(part,range)%i_part(ist))
                Col_FSC_seq  = mat(part,range)%i_seq(ist)
                gcspart = BCS_parts%icover_ibeg(mat(part,range)%i_part(ist))+Col_FSC_seq-1
                ! Debugging information
                if(iprint_DM>=5) write(*,5) myid,Col_FSC_part,Col_FSC_seq
                if(iprint_DM>=5) write(*,6) myid,BCS_parts%xcover(gcspart),&
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
                      ! Find processor that this block belongs to
                      proc = proc_block(SC_to_refx(SCblockr,SCblockc), &
                           SC_to_refy(SCblockr,SCblockc))
                      if(iprint_DM>=5) write(*,4) myid,proc,SCblockr,SCblockc
                      ! Create atom numbers and store posn in array at point
                      j = row-Distrib%firstrow(proc)+1  
                      k = (SCblockc-1)*block_size_c + SCrowc-firstcol(proc)+1
                      if(iprint_DM>=5) write(*,9) myid,proc,SCblockc,SCrowc,block_size_c,firstcol(proc)
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
!                      write(*,*) myid,' FSCpart, atom and xyz: ',FSCpart,parts%icell_beg(FSCpart)+mat(part,range)%i_seq(ist)-1,&
!                           x_atom_cell(parts%icell_beg(FSCpart)+mat(part,range)%i_seq(ist)-1),&
!                           y_atom_cell(parts%icell_beg(FSCpart)+mat(part,range)%i_seq(ist)-1),&
!                           z_atom_cell(parts%icell_beg(FSCpart)+mat(part,range)%i_seq(ist)-1)
!                      Distrib%loc(proc,j,k)%dx = BCS_parts%xcover(gcspart)- &
!                           x_atom_cell(parts%icell_beg(FSCpart)+mat(part,range)%i_seq(ist)-1)
!                      Distrib%loc(proc,j,k)%dy = BCS_parts%ycover(gcspart)- &
!                           y_atom_cell(parts%icell_beg(FSCpart)+mat(part,range)%i_seq(ist)-1)
!                      Distrib%loc(proc,j,k)%dz = BCS_parts%zcover(gcspart)- &
!                           z_atom_cell(parts%icell_beg(FSCpart)+mat(part,range)%i_seq(ist)-1)
!                      if(iprint_DM>=5) write(*,7) myid,proc,j,k,Distrib%loc(proc,j,k)
                      if(j>maxr.OR.k>maxc) write(*,*) 'Problem! j,k,max row,col: ',j,k,maxr,maxc
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
1   format(10x,'Processor: ',i5,' Partition: ',i5)
2   format(10x,'Processor: ',i5,' Atom: ',i5)
3   format(10x,'Processor: ',i5,' Neighbour: ',i5)
4   format(10x,'Processor: ',i5,' Recv: ',i5,' SC blocks: ',2i5)
5   format(10x,'Processor: ',i5,' FSC part, seq: ',2i5)
6   format(10x,'Processor: ',i5,' Neigh xyz: ',3f15.10)
7   format(10x,'Processor: ',i5,' Recv Proc: ',i5,' PS Atom j,k,loc: ',5i5,' dr: ',3f10.5)
8   format(10x,'Processor: ',i5,' Recv Proc: ',i5,' First row, col: ',2i5)
9   format(10x,'Processor: ',i5,' Recv Proc: ',i5,' SC block, atom, size: ',3i5,' firstcol: ',i5)
10  format(10x,'Processor: ',i5,' Maxr,c: ',2i5,' Maxrow,col: ',2i5)
11  format(10x,'Processor: ',i5,' PrepSend to : ',i5,' Atom, Supfn: ',2i5)
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
!!   datatypes, common, global_module, mpi, numbers, ScalapackFormat, maxima_module, GenComms
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
!!    Changed so that the matrix to be distributed is passed as an argument (allows more than
!!    one type of matrix to be sent) and added more comments
!!   16/04/2002 dave
!!    Added k-point to arguments passed
!!   16/04/2002 dave
!!    Added phase calculation to construction of Hamiltonian and made all appropriate variables complex
!!   18/04/2002 dave
!!    Added matrix we distribute into to list of arguments (i.e. allow distribution of more than one matrix)
!!    Also removed TODO entry, as I'd done it
!!   23/04/2002 dave
!!    Added intent to passed variables
!!   23/04/2002 dave
!!    Added DistributeData variable to carry details of matrix-specific arrays for distribution
!!   24/04/2002 dave
!!    Fixed a problem with accumulation of H for on-processor transfers
!!   08:18, 2003/07/28 dave
!!    Completely rewrote the accumulation loop to account for multiple CQ elements adding onto a single SC element.
!!    Changed Distrib structure, introduced loop over CQ elements, reworked various checks
!!   2004/11/10 drb
!!    Changed nsf to come from maxima, not common
!!   10:09, 13/02/2006 drb 
!!    Removed all explicit references to data_ variables and rewrote in terms of new 
!!    matrix routines
!!  SOURCE
!!
  subroutine DistributeCQ_to_SC(Distrib,matA,kk,SCmat)

    use datatypes
    use global_module, ONLY: numprocs, iprint_DM
    use mpi
    use numbers, ONLY: zero, minus_i
    use ScalapackFormat, ONLY: proc_start, block_size_r, block_size_c
    use GenComms, ONLY: my_barrier, myid
    use GenBlas, ONLY: copy
    use mult_module, ONLY: return_matrix_value_pos, matrix_pos
    
    implicit none

    ! Passed variables
    type(DistributeData), intent(in) :: Distrib
    integer, intent(in) :: matA
    real(double), intent(in) :: kk(3)
    complex(double_cplx), intent(out) :: SCmat(:,:)

    ! Local variables
    integer :: send_proc,recv_proc,send_size,recv_size
    integer :: sendtag, recvtag, stat
    integer :: srow_size,scol_size,rrow_size,rcol_size
    integer :: req, i, j, k, l, ierr, wheremat
    integer, dimension(MPI_STATUS_SIZE) :: mpi_stat
    real(double) :: phase,rfac,ifac,hreal,himag

    send_proc = myid
    recv_proc = myid
    ! Distribute data: loop over processors and issue sends and receives as appropriate
    SCmat = zero
    do i=1,numprocs
       if(iprint_DM>=4) write(*,1) myid,i,send_proc,recv_proc
       ! Sizes
       srow_size = Distrib%send_rows(send_proc+1)                    ! Sizes of remote 
       scol_size = proc_start(send_proc+1)%cols*block_size_c ! "chunk"
       rrow_size = Distrib%num_rows(recv_proc+1)                     ! Sizes of local
       rcol_size = proc_start(myid+1)%cols*block_size_c      ! "chunk"
       send_size = srow_size*scol_size
       recv_size = rrow_size*rcol_size
       ! Allocate memory
       allocate(SendBuffer(srow_size,scol_size),STAT=stat)
       if(stat/=0) call cq_abort("DiagModule: Can't alloc SendBuffer",stat)
       allocate(RecvBuffer(rrow_size,rcol_size),STAT=stat)
       if(stat/=0) call cq_abort("DiagModule: Can't alloc RecvBuffer",stat)
       ! Zero SendBuffer
       SendBuffer = cmplx(zero,zero,double_cplx)
       ! Zero RecvBuffer
       RecvBuffer = cmplx(zero,zero,double_cplx)
       ! On-site
       if(send_proc==myid.AND.recv_proc==myid) then 
          if(iprint_DM>=5) write(*,*) 'num_rows, send_rows: ',Distrib%num_rows(myid+1),Distrib%send_rows(myid+1)
          if(iprint_DM>=5) write(*,11) myid,srow_size,scol_size,send_size,rrow_size,rcol_size,recv_size
          ! Fill local copy of H matrix
          do j=1,srow_size
             do k=1,scol_size
                if(Distrib%images(send_proc+1,j,k)%n_elements>0) then
                   do l=1,Distrib%images(send_proc+1,j,k)%n_elements
                      wheremat = matrix_pos(matA,Distrib%images(send_proc+1,j,k)%where(l)%loci, &
                           Distrib%images(send_proc+1,j,k)%where(l)%locj,&
                           Distrib%images(send_proc+1,j,k)%where(l)%supfn_row,&
                           Distrib%images(send_proc+1,j,k)%where(l)%supfn_col)
                      if(iprint_DM>=5) write(*,7) myid,send_proc,j,k,wheremat
                      phase = kk(1)*Distrib%images(send_proc+1,j,k)%where(l)%dx + &
                           kk(2)*Distrib%images(send_proc+1,j,k)%where(l)%dy + &
                           kk(3)*Distrib%images(send_proc+1,j,k)%where(l)%dz
                      rfac = cos(phase)* return_matrix_value_pos(matA,wheremat)
                      ifac = sin(phase)* return_matrix_value_pos(matA,wheremat)
                      ! Care here - we need to accumulate
                      SCmat(Distrib%start_row(myid+1)+j-1,k) = SCmat(Distrib%start_row(myid+1)+j-1,k) + &
                           cmplx(rfac,ifac,double_cplx)
                   end do ! Distrib%images%n_elements
                end if ! n_elements>0
             end do ! k=1,scol_size
          end do ! j=1,srow_size
          if(iprint_DM>=5) write(*,8) myid,Distrib%start_row(recv_proc+1)
          if(iprint_DM>=4) write(*,*) '  Done on-proc'
          deallocate(RecvBuffer,STAT=stat)
          deallocate(SendBuffer,STAT=stat)
       ! Send and receive data to/from remote processors
       else ! if(send_proc==myid.AND.recv_proc==myid)
          ! ---------------
          ! Fill SendBuffer
          ! ---------------
          if(send_size>0) then
             do j=1,srow_size
                do k=1,scol_size
                   if(Distrib%images(send_proc+1,j,k)%n_elements>0) then
                      do l=1,Distrib%images(send_proc+1,j,k)%n_elements
                         wheremat = matrix_pos(matA,Distrib%images(send_proc+1,j,k)%where(l)%loci, &
                              Distrib%images(send_proc+1,j,k)%where(l)%locj,&
                              Distrib%images(send_proc+1,j,k)%where(l)%supfn_row,&
                              Distrib%images(send_proc+1,j,k)%where(l)%supfn_col)
                         if(iprint_DM>=5) write(*,7) myid,send_proc,j,k,wheremat
                         phase = kk(1)*Distrib%images(send_proc+1,j,k)%where(l)%dx + &
                              kk(2)*Distrib%images(send_proc+1,j,k)%where(l)%dy + &
                              kk(3)*Distrib%images(send_proc+1,j,k)%where(l)%dz
                         rfac = cos(phase)*return_matrix_value_pos(matA,wheremat)
                         ifac = sin(phase)*return_matrix_value_pos(matA,wheremat)
                         ! Accumulate the data
                         SendBuffer(j,k)= SendBuffer(j,k)+cmplx(rfac,ifac,double_cplx)
                      end do ! l=Distrib%images%n_elements
                   end if ! n_elements>0
                end do ! j=srow_size
             end do ! k=scol_size
          end if ! if send_size>0
          ! ---------------
          ! Do the transfer
          ! ---------------
          ! MPI tags
          sendtag = myid + send_proc*numprocs
          recvtag = recv_proc + myid*numprocs
          ! Debugging output
          if(iprint_DM>=5) then
             write(*,*) 'Proc: ',myid,' Sizes: ',send_size,recv_size
             write(*,10) i,myid,send_proc,recv_proc,srow_size, scol_size, rrow_size, rcol_size
             write(*,*) 'About to send...'
          end if
          ! Issue non-blocking send and then receive
          if(send_size>0) then
             call MPI_isend(SendBuffer,send_size,MPI_DOUBLE_COMPLEX,&
                  send_proc,sendtag,MPI_COMM_WORLD,req,ierr)
          end if
          if(recv_size>0) then
             call MPI_recv(RecvBuffer,recv_size,MPI_DOUBLE_COMPLEX,&
                  recv_proc,recvtag,MPI_COMM_WORLD,mpi_stat,ierr)
          end if
          ! More debugging output
          if(iprint_DM>=5) then
             write(*,8) myid,Distrib%start_row(recv_proc+1)
          end if
          if(iprint_DM>=4) write(*,2) myid,i
          ! Put the data in place from the receive buffer
          if(iprint_DM>=4) write(*,*) myid,'rowsize, start, col: ',rrow_size, Distrib%start_row(recv_proc+1), rcol_size
          if(rrow_size > 0) then
             do j=1,rrow_size
                SCmat(Distrib%start_row(recv_proc+1)+j-1,1:rcol_size) = RecvBuffer(j,1:rcol_size)
             end do
          end if ! (rrow_size > 0)
          ! Now wait for the non-blocking send to finish before deallocating !
          if(iprint_DM>=4) write(*,13) myid,i
          if(send_size>0) call MPI_Wait(req,mpi_stat,ierr)
          deallocate(RecvBuffer,STAT=stat)
          if(stat/=0) call cq_abort("DiagModule: Failed to dealloc buffer",stat)
          deallocate(SendBuffer,STAT=stat)
          if(stat/=0) call cq_abort("DiagModule: Failed to dealloc buffer",stat)
          if(iprint_DM>=4) write(*,12) myid,i
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
!!   DistributeSC_to_ref(Scalapack formatted eigenvector chunk, Conquest formatted eigenvectore chunk)
!!  PURPOSE
!!   Once a diagonalisation has been performed, we want eigenvector coefficients for every band for
!!   a given support function to be returned, in order of energy, to the processor responsible for 
!!   the support function.  In other words, we have \psi_n = \sum_{i\alpha} c^n_{i\alpha} \phi_{i\alpha}
!! (roughly speaking) and we want all n values of c^n_{i\alpha} on the processor responsible for i.
!!
!!   This is reasonably simple to accomplish - the data transfer is done just by reversing the calls of
!!   DistributeCQ_to_SC.  However, the data that arrives back is then ordered correctly in terms of rows
!!   but has its columns in SC format - we want them in reference format ! So we use mapy (I think) to do
!!   this.  For a chunk which has come back, we know that we have all columns, so we work out which row
!!   block we're in (for a given row) and loop over columns, mapping to reference block and hence position
!!   in chunk.
!!  INPUTS
!!   real(double), dimension(:,:) :: SCeig - the piece of eigenvector matrix created by Scalapack
!!   real(double), dimension(:,:) :: localEig - local storage for eigenvector coefficients for atoms
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   17/04/2002
!!  MODIFICATION HISTORY
!!   18/04/2002 dave
!!    More documentation and finished off decoding of RecvBuffer (I think)
!!   23/04/2002 dave
!!    Added intents
!!   23/04/2002 dave
!!    Added DistributeData derived type to carry information about matrix-specific data - though actually here
!!    I think that the only bits used (firstrow) are matrix-independent.
!!   2004/11/10 drb
!!    Changed nsf to come from maxima, not common
!!   2006/11/02 16:59 dave
!!    Bug fix: only call MPI_Wait if MPI_isend called (i.e. send_size greater than zero)
!!  SOURCE
!!
  subroutine DistributeSC_to_ref(Distrib,SCeig,localEig)
  
    use datatypes
    use global_module, ONLY: numprocs, iprint_DM
    use mpi
    use numbers, ONLY: zero, minus_i
    use ScalapackFormat, ONLY: proc_start, block_size_r, block_size_c, mapy
    use GenComms, ONLY: my_barrier, myid

    implicit none

    ! Passed variables
    type(DistributeData), intent(in) :: Distrib
    complex(double_cplx), dimension(:,:), intent(in)  :: SCeig
    complex(double_cplx), dimension(:,:), intent(out) :: localEig

    ! Local variables
    integer :: send_proc,recv_proc,send_size,recv_size
    integer :: sendtag, recvtag, stat,rblock,cblock,refblock,roff,coff,req1,req2,ierr
    integer :: srow_size,scol_size,rrow_size,rcol_size
    integer, dimension(MPI_STATUS_SIZE) :: mpi_stat
    integer :: i,j,k

    send_proc = myid
    recv_proc = myid
    ! Distribute data: loop over processors and issue sends and receives as appropriate
    localEig = zero
    do i=1,numprocs
       if(iprint_DM>=4) write(*,1) myid,i,send_proc,recv_proc
       ! Sizes
       rrow_size = Distrib%send_rows(recv_proc+1)                    ! Sizes of remote 
       rcol_size = proc_start(recv_proc+1)%cols*block_size_c         ! "chunk"
       srow_size = Distrib%num_rows(send_proc+1)                     ! Sizes of local
       scol_size = proc_start(myid+1)%cols*block_size_c              ! "chunk"
       send_size = srow_size*scol_size
       recv_size = rrow_size*rcol_size
       ! Allocate memory
       allocate(SendBuffer(srow_size,scol_size),STAT=stat)
       if(stat/=0) call cq_abort("DiagModule: Can't alloc SendBuffer",stat)
       allocate(RecvBuffer(rrow_size,rcol_size),STAT=stat)
       if(stat/=0) call cq_abort("DiagModule: Can't alloc RecvBuffer",stat)
       ! Zero SendBuffer
       SendBuffer = cmplx(zero,zero,double_cplx)
       ! Zero RecvBuffer
       RecvBuffer = cmplx(zero,zero,double_cplx)
       if(iprint_DM>=5) write(*,8) myid,Distrib%start_row(recv_proc+1)
       ! On-site
       if(send_proc==myid.AND.recv_proc==myid) then 
          if(iprint_DM>=5) write(*,*) 'num_rows, send_rows: ',Distrib%num_rows(myid+1),Distrib%send_rows(myid+1)
          if(iprint_DM>=5) write(*,11) myid,srow_size,scol_size,send_size,rrow_size,rcol_size,recv_size
          ! Fill send buffer
          ! Normally we'd then send this, and copy our receive buffer out
          do j=1,srow_size
!             SendBuffer(j,1:scol_size) = SCeig(start_row(send_proc+1)+j-1,1:scol_size)
             RecvBuffer(j,1:scol_size) = SCeig(Distrib%start_row(send_proc+1)+j-1,1:scol_size)
          end do
          roff = Distrib%start_row(send_proc+1)
          ! Decode receive buffer into local eigenvector store
          do j=1,rrow_size
             rblock = aint(real((roff-1+j-1)/block_size_r))+1
             do k=1,rcol_size,block_size_c
                cblock = aint(real((k-1)/block_size_c))+1
                refblock = mapy(recv_proc+1,rblock,cblock)
                coff = (refblock-1)*block_size_c + 1
                if(iprint_DM>=5) write(*,3) myid,j,k,rblock,cblock,refblock,coff,Distrib%firstrow(recv_proc+1),RecvBuffer(j,k)
!                localEig(Distrib%firstrow(recv_proc+1)+j-1,coff:coff+block_size_c-1) = RecvBuffer(j,k:k+block_size_c-1)
                localEig(coff:coff+block_size_c-1,Distrib%firstrow(recv_proc+1)+j-1) = RecvBuffer(j,k:k+block_size_c-1)
             end do
          end do
          if(iprint_DM>=4) write(*,*) '  Done on-proc'
          deallocate(RecvBuffer,STAT=stat)
          deallocate(SendBuffer,STAT=stat)
       ! Send and receive data to/from remote processors
       else ! if(send_proc==myid.AND.recv_proc==myid)
          ! ---------------
          ! Fill SendBuffer
          ! ---------------
          if(send_size>0) then
             do j=1,srow_size
                SendBuffer(j,1:scol_size) = SCeig(Distrib%start_row(send_proc+1)+j-1,1:scol_size)
             end do
          end if ! if send_size>0
          ! ---------------
          ! Do the transfer
          ! ---------------
          ! MPI tags
          sendtag = myid + send_proc*2*numprocs
          recvtag = recv_proc + myid*2*numprocs
          ! Debugging output
          if(iprint_DM>=5) then
             write(*,*) 'Proc: ',myid,' Sizes: ',send_size,recv_size
             write(*,10) i,myid,send_proc,recv_proc,srow_size, scol_size, rrow_size, rcol_size
             write(*,*) 'About to send...'
          end if
          ! Issue non-blocking send and then receive
          if(send_size>0) then
             call MPI_isend(Distrib%start_row(send_proc+1),1,MPI_INTEGER,send_proc,sendtag,MPI_COMM_WORLD,req1,ierr)
             call MPI_isend(SendBuffer,send_size,MPI_DOUBLE_COMPLEX,send_proc,sendtag+1,MPI_COMM_WORLD,req2,ierr)
          end if
          if(recv_size>0) then
             call MPI_recv(roff,1,MPI_INTEGER,recv_proc,recvtag,MPI_COMM_WORLD,mpi_stat,ierr)
             call MPI_recv(RecvBuffer,recv_size,MPI_DOUBLE_COMPLEX,&
                  recv_proc,recvtag+1,MPI_COMM_WORLD,mpi_stat,ierr)
          end if
          if(iprint_DM>=4) write(*,2) myid,roff
          ! Put the data in place from the receive buffer
          if(rrow_size > 0) then
             do j=1,rrow_size
                rblock = aint(real((roff-1+j-1)/block_size_r))+1
                do k=1,rcol_size,block_size_c
                   cblock = aint(real((k-1)/block_size_c))+1
                   refblock = mapy(recv_proc+1,rblock,cblock)
                   coff = (refblock-1)*block_size_c + 1
                   if(iprint_DM>=5) write(*,3) myid,j,k,rblock,cblock,refblock,coff,Distrib%firstrow(recv_proc+1),RecvBuffer(j,k)
                   !                localEig(Distrib%firstrow(recv_proc+1)+j-1,coff:coff+block_size_c-1) = RecvBuffer(j,k:k+block_size_c-1)
                   localEig(coff:coff+block_size_c-1,Distrib%firstrow(recv_proc+1)+j-1) = RecvBuffer(j,k:k+block_size_c-1)
                end do
             end do
          end if ! (rrow_size > 0)
          ! Now wait for the non-blocking send to finish before deallocating !
          if(send_size>0) then
             call MPI_Wait(req1,mpi_stat,ierr)
             call MPI_Wait(req2,mpi_stat,ierr)
          end if
          deallocate(RecvBuffer,STAT=stat)
          if(stat/=0) call cq_abort("DiagModule: Failed to dealloc buffer",stat)
          deallocate(SendBuffer,STAT=stat)
          if(stat/=0) call cq_abort("DiagModule: Failed to dealloc buffer",stat)
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

!!****f* DiagModule/findFermi *
!!
!!  NAME 
!!   findFermi
!!  USAGE
!!   findFermi(eig,nbands,nkp,Ef)
!!   findFermi(eigenstates, number of bands, number of k points, fermi energy)
!!  PURPOSE
!!   Finds the fermi level given a set of eigenvalues at a number of k points
!!  INPUTS
!!   integer :: nbands, nkp - Numbers of bands and k points
!!   real(double), dimension(nbands,nkp), intent(in) :: eig - Eigenvalues for each k point
!!   real(double), intent(out) :: Ef - fermi energy
!!  USES
!!   datatypes, numbers, common
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   23/04/2002
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine findFermi(electrons,eig,nbands,nkp,Ef)

    use datatypes
    use numbers, ONLY : half,one
    use global_module, ONLY: iprint_DM
    use GenComms, ONLY: myid
    
    implicit none

    ! Passed variables
    integer :: nbands, nkp
    real(double), dimension(nbands,nkp), intent(in) :: eig
    real(double), intent(out) :: Ef
    real(double), intent(in) :: electrons

    ! Local variables
    real(double) :: thisElec, lowElec, highElec
    real(double) :: lowEf, highEf
    real(double), parameter :: incEf = one
    real(double), parameter :: tolElec = 1.0e-6_double

    if(iprint_DM>=2) write(*,5) myid,electrons
    ! Take first guess as double filling each band at first k point
    Ef = eig(int(electrons/2),1)
    call occupy(occ,eig,Ef,thisElec,nbands,nkp)
    ! Find two values than bracket true Ef
    if(thisElec<electrons) then ! We've found a lower bound
       if(iprint_DM>=4) write(*,3) myid,Ef
       lowEf = Ef
       lowElec = thisElec
       highEf = lowEf + incEf
       call occupy(occ,eig,highEf,highElec,nbands,nkp)
       do while(highElec<electrons) ! Increase upper bound
          lowEf = highEf
          lowElec = highElec
          highEf = highEf + incEf
          call occupy(occ,eig,highEf,highElec,nbands,nkp)
          if(iprint_DM>=4) write(*,6) myid,highEf,highElec
       end do
    else ! We have an upper bound
       if(iprint_DM>=4) write(*,4) myid,Ef
       highEf = Ef
       highElec = thisElec
       lowEf = highEf - incEf
       call occupy(occ,eig,lowEf,lowElec,nbands,nkp)
       do while(lowElec>electrons) ! Decrease lower bound
          highEf = lowEf
          highElec = lowElec
          lowEf = lowEf - incEf
          call occupy(occ,eig,lowEf,lowElec,nbands,nkp)
          if(iprint_DM>=4) write(*,6) myid,lowEf,lowElec
       end do
    end if
    if(iprint_DM>=3) write(*,2) myid,lowEf,highEf
    Ef = half*(lowEf + highEf)
    call occupy(occ,eig,Ef,thisElec,nbands,nkp)
    do while(abs(thisElec - electrons)>tolElec)
       if(thisElec>electrons) then
          highElec = thisElec
          highEf = Ef
       else
          lowElec = thisElec
          lowEf = Ef
       end if
       Ef = half*(lowEf + highEf)
       call occupy(occ,eig,Ef,thisElec,nbands,nkp)
    end do
    if(iprint_DM>=2) write(*,1) Ef
    return
1   format(10x,'Fermi level is ',f12.5)
2   format(10x,'Proc: ',i5,' bracketed Ef: ',2f12.5)
3   format(10x,'Proc: ',i5,' In findFermi, found lower bound',f12.5)
4   format(10x,'Proc: ',i5,' In findFermi, found upper bound',f12.5)
5   format(10x,'Proc: ',i5,' In findFermi, searching for Ne: ',f12.5)
6   format(10x,'Proc: ',i5,' In findFermi, level, Ne: ',2f12.5)
  end subroutine findFermi
!!***

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
!!   Populates the eigenstates up to a given fermi level at each k-point.  Called by findFermi to 
!!   find the fermi level.
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
!!
!!  SOURCE
!!
  subroutine occupy(occ,ebands,Ef,electrons,nbands,nkp)

    use datatypes
    use numbers, ONLY: zero, two
    use global_module, ONLY: iprint_DM
    use GenComms, ONLY: myid

    implicit none

    ! Passed variables
    integer, intent(in) :: nbands, nkp
    real(double), dimension(nbands,nkp), intent(out) :: occ
    real(double), dimension(nbands,nkp), intent(in) :: ebands
    real(double), intent(in) :: Ef
    real(double), intent(out) :: electrons

    ! Local variables
    integer :: ikp, iband

    electrons = zero
    do ikp = 1,nkp
       do iband = 1,nbands
          occ(iband,ikp) = two*wtk(ikp)*fermi(ebands(iband,ikp)-Ef,kT)
          electrons = electrons + occ(iband,ikp)
       end do
    end do
    if(iprint_DM>=5) write(*,1) myid,Ef,electrons
    return
1   format(10x,'In occupy on proc: ',i5,' For Ef of ',f8.5,' we get ',f12.5,' electrons')
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
!!    Small fix to prevent maths overflows by only calculating exponential if x well bounded
!!  SOURCE
!!
  real(double) function fermi(E,kT)

    use datatypes
    use numbers, ONLY: zero, one, half
    
    implicit none

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

!!****f* DiagModule/buildK *
!!
!!  NAME 
!!   buildK - makes K 
!!  USAGE
!! 
!!  PURPOSE
!!   Builds K from eigenvectors - this involves working out which processors we're going to need 
!!   which eigenvectors from, fetching the data and building the matrix
!!
!!   N.B. The conjugation of one set of eigenvector coefficients takes place when the dot product is
!!   performed - dot maps onto zdotc which conjugates the FIRST vector.  Really we should conjugate the
!!   second, but as K is real, it shouldn't matter.
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
!!    This routine is now working, so tidied - added deallocate statements, iprint_DM levels to write
!!    statements and deleted unnecessary rubbish
!!   2004/11/10 drb
!!    Changed nsf to come from maxima, not common
!!   10:09, 13/02/2006 drb 
!!    Removed all explicit references to data_ variables and rewrote in terms of new 
!!    matrix routines
!!   2006/10/02 17:52 dave
!!    Added deallocate for norb_send (thanks TM, TO)
!!  SOURCE
!!
  subroutine buildK(range, matA, occs, kps, weight, localEig)

    !use maxima_module, ONLY: mx_nponn, mx_at_prim
    use matrix_module, ONLY: matrix, matrix_halo
    use group_module, ONLY: parts
    use primary_module, ONLY: bundle
    use cover_module, ONLY: BCS_parts
    use ScalapackFormat, ONLY : CC_to_SC,maxrow,maxcol,proc_block, &
         SC_to_refx,SC_to_refy, block_size_r, block_size_c, blocks_c, proc_start, matrix_size
    use global_module, ONLY: numprocs, iprint_DM, id_glob, ni_in_cell, x_atom_cell, y_atom_cell, z_atom_cell
    use numbers, ONLY: very_small
    use mpi
    use GenBlas, ONLY: dot
    use GenComms, ONLY: myid
    use mult_module, ONLY: store_matrix_value_pos, matrix_pos
    use matrix_data, ONLY: mat, halo
    use species_module, ONLY: nsf_species

    implicit none

    ! Passed variables
    real(double), dimension(:), intent(in) :: occs
    real(double), dimension(3), intent(in) :: kps
    real(double) :: weight
    integer :: matA, range
    complex(double_cplx), dimension(:,:), intent(in) :: localEig
    ! Local variables
    type(Krecv_data), dimension(:), allocatable :: recv_info
    integer :: part, memb, neigh, ist, prim_atom, owning_proc, locatom
    integer :: Row_FSC_part, Row_FSC_seq, Col_FSC_part, Col_FSC_seq, FSC_atom
    integer :: SCblockr, SCblockc, SCrowc, i, j, k, proc, stat, supfn_r, supfn_c
    integer :: maxloc, maxint, maxsend, curr, gcspart, CC, orb_count
    integer :: len, send_size, recv_size, send_proc, recv_proc, nsf1, sendtag, recvtag
    integer :: req1, req2, ierr, atom, inter, prim, wheremat, row_sup, col_sup
    integer, dimension(:,:), allocatable :: ints, atom_list, send_prim, send_info, send_orbs, send_off
    integer, dimension(:), allocatable :: current_loc_atoms, LocalAtom, num_send, norb_send, &
         send_FSC, recv_to_FSC, mapchunk, prim_orbs
    integer, dimension(MPI_STATUS_SIZE) :: mpi_stat
    real(double) :: phase, rfac, ifac, rcc, icc, rsum
    complex(double_cplx) :: zsum
    complex(double_cplx), dimension(:,:), allocatable :: RecvBuffer, SendBuffer
    logical :: flag
    integer :: FSCpart, ipart

    if(iprint_DM>=2.AND.myid==0) write(*,fmt='(10x,"Entering buildK ",i4)') matA
    ! Allocate data and zero arrays
    allocate(ints(numprocs,bundle%mx_iprim),current_loc_atoms(numprocs),atom_list(numprocs,bundle%mx_iprim),&
         LocalAtom(ni_in_cell),send_prim(numprocs,bundle%mx_iprim),&
         num_send(numprocs),norb_send(numprocs),STAT=stat)
    if(stat/=0) call cq_abort('buildK: Error allocating ints, current_loc_atoms and atom_list !',stat)
    ints = 0
    current_loc_atoms = 0
    atom_list = 0
    LocalAtom = 0
    send_prim = 0
    num_send = 0
    norb_send = 0
    if(iprint_DM>=3) write(*,*) 'buildK: Stage one'
!    write(*,*) myid,' send_prim: ',send_prim
    ! Step one - work out which processors we need to exchange data with
    do part = 1,bundle%groups_on_node ! Loop over primary set partitions
       if(iprint_DM>=5) write(*,1) myid,part
       if(bundle%nm_nodgroup(part)>0) then ! If there are atoms in partition
          CC = parts%ngnode(parts%inode_beg(myid+1)+part-1)
          do memb = 1,bundle%nm_nodgroup(part) ! Loop over atoms
             if(iprint_DM>=5) write(*,2) myid,memb
             prim_atom = bundle%nm_nodbeg(part)+memb-1
             do neigh = 1, mat(part,range)%n_nab(memb) ! Loop over neighbours of atom
                if(iprint_DM>=5) write(*,3) myid,neigh
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
                if(iprint_DM>=5) write(*,*) 'prim, neigh, FSC: ',prim_atom, neigh, FSC_atom
                if(iprint_DM>=5) write(*,*) 'curr_loc_atoms: ',current_loc_atoms(owning_proc)
                if(current_loc_atoms(owning_proc)>0) then
                   do i=1,current_loc_atoms(owning_proc)
                      if(atom_list(owning_proc,i)==FSC_atom) then
                         if(iprint_DM>=5) write(*,*) 'Loc atom: ',i, LocalAtom(FSC_atom)
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
    if(iprint_DM>=3) write(*,*) 'buildK: Stage two'
    maxloc = 0
    maxint = 0
    maxsend = 0
    do i=1,numprocs
       if(iprint_DM>=4) write(*,*) myid,' Curr loc atoms: ',i,current_loc_atoms(i)
       if(current_loc_atoms(i)>maxloc) maxloc = current_loc_atoms(i)
       do j=1,bundle%mx_iprim ! Needs to be mx_iprim because goes over primary atoms on REMOTE processors
          if(ints(i,j)>maxint) maxint = ints(i,j)
          if(send_prim(i,j)>0) num_send(i) = num_send(i) + 1
          norb_send(i) = norb_send(i) + send_prim(i,j)
          if(iprint_DM>=5) write(*,4) myid,j,send_prim(i,j),num_send(i)
       end do
       if(num_send(i)>maxsend) maxsend = num_send(i)
    end do
    if(iprint_DM>=4) write(*,*) myid,' Maxima: ',maxloc, maxint, maxsend
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
       if(iprint_DM>=5) write(*,1) myid,part
       if(bundle%nm_nodgroup(part)>0) then ! If there are atoms in partition
          CC = parts%ngnode(parts%inode_beg(myid+1)+part-1)
          do memb = 1,bundle%nm_nodgroup(part) ! Loop over atoms
             if(iprint_DM>=5) write(*,2) myid,memb
             prim_atom = bundle%nm_nodbeg(part)+memb-1
             do neigh = 1, mat(part,range)%n_nab(memb) ! Loop over neighbours of atom
                if(iprint_DM>=5) write(*,3) myid,neigh
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
                if(iprint_DM>=5) write(*,*) myid,' own, FSC, loc: ',owning_proc, FSC_atom, locatom, &
                     recv_info(owning_proc)%ints(locatom)
                recv_info(owning_proc)%ints(locatom) = recv_info(owning_proc)%ints(locatom) + 1
                if(iprint_DM>=5) write(*,*) myid,' ints: ',recv_info(owning_proc)%ints(locatom)
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
    len = 0
    do i=1,matrix_size
       if(abs(occs(i))>very_small) then
          len = len+1
          if(myid==0.AND.iprint_DM>=4) write(*,*) 'Occ is ',occs(i)
       end if
    end do
    if(iprint_DM>=3) write(*,*) 'buildK: Stage three len:',len, matA
    ! Step three - loop over processors, send and recv data and build K
    allocate(send_fsc(bundle%mx_iprim),recv_to_FSC(bundle%mx_iprim),mapchunk(bundle%mx_iprim),STAT=stat)
    if(stat/=0) call cq_abort('buildK: Error allocating send_fsc, recv_to_FSC and mapchunk',stat)
    send_fsc = 0
    recv_to_FSC = 0
    mapchunk = 0
    send_proc = myid
    recv_proc = myid
    do i=1,numprocs
       send_size = len*norb_send(send_proc+1)!num_send(send_proc+1)*nsf
       recv_size = len*recv_info(recv_proc+1)%orbs!current_loc_atoms(recv_proc+1)*nsf
       if(iprint_DM>=4) write(*,*) 'Send and recv sizes: ',send_size, recv_size
       ! Fill SendBuffer
       allocate(SendBuffer(len,norb_send(send_proc+1)),STAT=stat)
       if(stat/=0) call cq_abort('buildK: Unable to allocate SendBuffer !',stat)
       if(iprint_DM>=4) write(*,*) 'Filling SendBuffer'
       orb_count = 0
       do j=1,num_send(send_proc+1)
          do nsf1=1,send_orbs(send_proc+1,j)
             orb_count = orb_count+1
             SendBuffer(1:len,orb_count) = localEig(1:len,send_off(send_proc+1,j)+nsf1)
          end do
          ! We also need to send a list of what FSC each primary atom sent corresponds to - use bundle%ig_prim
         send_FSC(j) = bundle%ig_prim(send_info(send_proc+1,j))
         if(iprint_DM>=4) write(*,*) 'Building send_FSC: ',send_info(send_proc+1,j), &
              bundle%ig_prim(send_info(send_proc+1,j)),send_FSC(j)
       end do
       if(orb_count/=norb_send(send_proc+1)) call cq_abort("Orbital mismatch in buildK: ",orb_count,norb_send(send_proc+1))
       sendtag = myid + send_proc*2*numprocs
       recvtag = recv_proc + myid*2*numprocs
       if(iprint_DM>=4) write(*,*) 'Sending'
       ! Now send
       if(send_size>0) then
          if(send_proc/=myid) then
             call MPI_issend(send_FSC,num_send(send_proc+1),MPI_INTEGER,send_proc,sendtag,MPI_COMM_WORLD,req1,ierr)
             call MPI_issend(SendBuffer,send_size,MPI_DOUBLE_COMPLEX,send_proc,sendtag+1,MPI_COMM_WORLD,req2,ierr)
          end if
       end if
       ! Now receive data
       if(iprint_DM>=4) write(*,*) 'Alloc RecvBuffer ',len,recv_info(recv_proc+1)%orbs
       !allocate(RecvBuffer(len,current_loc_atoms(recv_proc+1)*nsf),STAT=stat)
       allocate(RecvBuffer(len,recv_info(recv_proc+1)%orbs),STAT=stat)
       if(stat/=0) call cq_abort('buildK: Unable to allocate RecvBuffer !',stat)
       if(iprint_DM>=4) write(*,*) 'Recving'
       if(recv_size>0) then
          if(recv_proc/=myid) then
             call MPI_recv(recv_to_FSC,current_loc_atoms(recv_proc+1),MPI_INTEGER,recv_proc,recvtag,MPI_COMM_WORLD,mpi_stat,ierr)
             if(iprint_DM>=4) write(*,*) 'Got recv_to_FSC'
             call MPI_recv(RecvBuffer,recv_size,MPI_DOUBLE_COMPLEX,&
                  recv_proc,recvtag+1,MPI_COMM_WORLD,mpi_stat,ierr)
             if(iprint_DM>=4) write(*,*) 'Got RecvBuffer'
          else
             if(iprint_DM>=4) write(*,*) 'On-proc: getting recv_to_FSC'
             recv_to_FSC(1:current_loc_atoms(recv_proc+1)) = send_FSC(1:current_loc_atoms(recv_proc+1))
             if(iprint_DM>=4) write(*,*) 'On-proc: getting RecvBuffer'
             RecvBuffer(1:len,1:recv_info(recv_proc+1)%orbs) = SendBuffer(1:len,1:recv_info(recv_proc+1)%orbs)
          end if
          if(iprint_DM>=4) write(*,*) 'Doing the mapchunk', recv_to_FSC
          do j=1,current_loc_atoms(recv_proc+1)
             mapchunk(j) = LocalAtom(recv_to_FSC(j))
          end do
          if(iprint_DM>=4) write(*,*) 'filling buffer'
          do j=1,len
             RecvBuffer(j,1:recv_info(recv_proc+1)%orbs) = RecvBuffer(j,1:recv_info(recv_proc+1)%orbs)*0.5_double*occs(j)
          end do
          orb_count = 0
          do atom = 1,current_loc_atoms(recv_proc+1)
             locatom = mapchunk(atom)
             if(iprint_DM>=4) write(*,*) 'Atom, loc: ',atom,locatom,recv_info(recv_proc+1)%ints(locatom)
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
                if(iprint_DM>=4) write(*,*) 'Inter: ',inter,prim
                phase = kps(1)*recv_info(recv_proc+1)%dx(inter,locatom) + kps(2)*recv_info(recv_proc+1)%dy(inter,locatom) + &
                     kps(3)*recv_info(recv_proc+1)%dz(inter,locatom)
                if(iprint_DM>=5) write(*,*) 'Prim, where, phase: ',prim, whereMat, phase
                rfac = cos(phase)
                ifac = sin(phase)
                do row_sup = 1,recv_info(recv_proc+1)%ndimj(locatom)
                   do col_sup = 1,nsf_species(bundle%species(prim))
                      whereMat = matrix_pos(matA,recv_info(recv_proc+1)%prim_atom(inter,locatom), &
                           recv_info(recv_proc+1)%locj(inter,locatom),col_sup,row_sup) 
                      zsum = dot(len,localEig(1:len,prim_orbs(prim)+col_sup),1,RecvBuffer(1:len,orb_count+row_sup),1)
                      call store_matrix_value_pos(matA,whereMat,real(zsum*cmplx(rfac,ifac,double_cplx),double))
                   end do ! col_sup=nsf
                end do ! row_sup=nsf
             end do ! inter=recv_info%ints
             ! Careful - we only want to increment after ALL interactions done
             orb_count = orb_count + recv_info(recv_proc+1)%ndimj(locatom)
          end do ! atom=current_loc_atoms
       end if ! recv_size>0
       if(iprint_DM>=4) write(*,*) 'Calling MPI_Wait'
       if(send_size>0.AND.myid/=send_proc) then
          call MPI_Wait(req1,mpi_stat,ierr)
          call MPI_Wait(req2,mpi_stat,ierr)
       end if
       if(iprint_DM>=4) write(*,*) 'Calling dealloc'
       deallocate(RecvBuffer,STAT=stat)
       if(stat/=0) call cq_abort("buildK: Failed to dealloc buffer",stat) 
       if(iprint_DM>=4) write(*,*) 'Calling dealloc'
       deallocate(SendBuffer,STAT=stat)
       if(stat/=0) call cq_abort("buildK: Failed to dealloc buffer",stat)
       ! Increment/decrement recv and send, and wrap
       ! Remember that we go from 0->numprocs-1
       if(iprint_DM>=4) write(*,*) 'Doing proc thang'
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
    !    write(*,*) 'Returning'
    return
1   format(10x,'Processor: ',i5,' Partition: ',i5)
2   format(10x,'Processor: ',i5,' Atom: ',i5)
3   format(10x,'Processor: ',i5,' Neighbour: ',i5)
4   format(10x,'Proc: ',i5,' Prim, send_prim, num_send: ',3i5)
  end subroutine buildK
!!***

end module DiagModule
