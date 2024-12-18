!-*- mode: F90; mode: font-lock -*-
! -----------------------------------------------------------------------------
! $Id$
! -----------------------------------------------------------------------------
! Module ScalapackFormat
! -----------------------------------------------------------------------------
! Code area 4: density matrix
! -----------------------------------------------------------------------------

!!****h* Conquest/ScalapackFormat *
!!  NAME
!!   ScalapackFormat
!!  PURPOSE
!!   Within Conquest and ParaDens, we want to introduce parallel   
!!   diagonalisation using ScaLAPACK.  This small module contains
!!   subroutines to create the different formats required: reference 
!! (in other words the "unpacked" or global ScaLAPACK matrix), Scalapack 
!! (the matrix as distributed over processors for ScaLAPACK)   
!!   and the relationship to Conquest primary sets and atom      
!!   numbers.                                                    
!!
!!   Scalapack has two basic formats: reference and Scalapack. We
!!   need to know (a) how these relate to each other and (b) how they
!!   relate to the matrices as stored by Conquest.  The reference 
!!   format is nominal - it's how the matrix would look on a single
!!   processor machine (i.e. written out globally).  The Scalapack is
!!   necessary - it's how the data is actually distributed when the 
!!   work is done on it by Scalapack.  As explained in the notes for
!!   Conquest diagonalisation, Scalapack uses the "2D block cyclic"
!!   distribution.  We insist that when the matrix is written in 
!!   Scalapack form (i.e. as a big matrix but with blocks grouped by
!!   processor) the rows are ordered exactly as in Conquest processor/
!!   primary set order - so that we only have to reorganise by column
!!   and then send.
!!  USES
!!   common, GenComms, global_module, group_module, maxima_module
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   28th January 2002
!!  MODIFICATION HISTORY
!!   07/02/2002 dave (Happy Birthday Robert)
!!    Tidied, renamed module, added documentation and Robodoc headers
!!   08/02/2002 dave
!!    Changed names of conversion arrays (ref_to_SC and SC_to_ref)
!!   11/02/2002 dave
!!    Removed extraneous print statements
!!   19/02/2002 dave (Happy Birthday Christopher)
!!    Imported to ParaDens
!!   22/02/2002 dave
!!    Added processor_data type (may need rethinking later) to store
!!    the starting point for the data for each processor in SC format
!!   28/02/2002 dave
!!    Added an array to convert from CC partition and sequence to SC column
!!   19/03/2002 dave
!!    Added CQ_to_SC which stores all Conquest indexing for an atom equivalent
!!    to a matrix row; also improved the arrangement of many routines and added
!!    support function numbers to each of them
!!   05/04/2002 dave
!!    Many further changes mainly related to the change from atom-by-atom indexing
!!    to row-by-row indexing (which is much more general and is really a requirement)
!!   16/04/2002 dave
!!    Moved output about start of subroutine from DiagModule to individual subroutines
!!   17/06/2002 dave
!!    Improved headers and tidied
!!   16:58, 10/11/2004 dave 
!!    Removed inappropriate common nsf declaration
!!   2008/02/04 08:23 dave
!!    Changed for output to file not stdout
!!   2008/05/19 ast
!!    Added timers
!!   2011/02/13 L.Tong
!!    Added k-point parallelisation
!!   2014/09/15 18:30 lat
!!    fixed call start/stop_timer to timer_module (not timer_stdlocks_module !)
!!   2017/06/22 13:47 dave
!!    Changes alongside DiagModule to split initialisation
!!  SOURCE
!!
module ScalapackFormat

  use global_module,          only: io_lun, area_DM, iprint_DM
  use timer_module,           only: start_timer, stop_timer
  use timer_stdclocks_module, only: tmr_std_allocation

  implicit none

  save
!!***

!!****s* ScalapackFormat/block_atom *
!!  NAME
!!   block_atom
!!  PURPOSE
!!   Stores data for atoms in Scalapack/reference/Conquest
!!    transformation arrays
!!  AUTHOR
!!   D.R.Bowler
!!  SOURCE
!!
  type block_atom
     integer :: part, atom, support_fn
  end type block_atom
!!***

!!****s* ScalapackFormat/processor_data *
!!  NAME
!!   processor_data
!!  PURPOSE
!!   Contains data about processors and where their data lies - 
!!    in particular where in the Scalapack rows their atoms start
!!  AUTHOR
!!   D.R.Bowler
!!  SOURCE
!!
  type processor_data
     integer :: block, atom
     integer :: rows, cols
     integer :: startrow
  end type processor_data
!!***

!!****s* ScalapackFormat/CC_to_SC_map *
!!  NAME
!!   CC_to_SC_map
!!  PURPOSE
!!   Maps from the fundamental simulation cell partition (CC)/sequence
!!   labelling to the Scalapack column block/atom labelling
!!  AUTHOR
!!   D.R.Bowler
!!  SOURCE
!!
  type CC_to_SC_map
     integer :: block_c, row_c
     integer :: block_r, row_r
  end type CC_to_SC_map
!!***

!!****s* ScalapackFormat/CQ_to_SC_map *
!!  NAME
!!   CQ_to_SC_map
!!  PURPOSE
!!   Gives for each row in matrix the Conquest indexing labels
!!  AUTHOR
!!   D.R.Bowler
!!  SOURCE
!!
  type CQ_to_SC_map
     integer :: support_fn
     integer :: atom
     integer :: partition
     integer :: CClabel
  end type CQ_to_SC_map
!!***

  ! Dimension of matrix (always square)
  integer :: matrix_size
  ! Dimension of padded H (and S) matrix 
  logical :: flag_padH
  integer :: matrix_size_padH
  ! Block sizes (not necessarily square)
  integer :: block_size_r, block_size_c
  ! Processor group
  integer :: proc_groups, pgid  ! proc_goups = tota number of process-groups and pgid = the local group id
  integer, dimension(:), allocatable :: pgroup ! proc_group index of a given process
  integer, dimension(:), allocatable :: N_procs_in_pg  ! number of processes in a given group
  integer, dimension(:), allocatable :: N_kpoints_in_pg  ! number of k-points in a given group
  integer, dimension(:,:), allocatable :: pg_procs  ! given a proc_group, stores the list of proc ids in the group
  integer, dimension(:,:), allocatable :: pg_kpoints  ! given a proc_group, strores the list of k-points in the group
  integer :: nprocs_max, nkpoints_max  ! max number of procs and kpoints can be in a group
  ! Processor grid size
  integer :: proc_rows, proc_cols
  ! Numbers of blocks, max blocks for processor, processors
  integer :: blocks_r,blocks_c,maxrow,maxcol

  ! Processor data
  type(processor_data), dimension(:), allocatable :: proc_start
  ! FSC to SC map
  type(CC_to_SC_map), dimension(:,:,:), allocatable :: CC_to_SC
  type(CQ_to_SC_map), dimension(:), allocatable :: CQ2SC_row_info
  ! Mask for my rows in SC format
  integer, dimension(:), allocatable :: my_row

  ! -------------------------------------------------------
  ! Allocatable variables
  ! -------------------------------------------------------
  ! Map from processor grid to processors, for each k-point process-group
  integer, allocatable, dimension(:,:,:) :: procid
  ! On-processor map giving reference position for a block
  ! (i.e. for any block on a processor, give x and y in ref format)
  integer, allocatable, dimension(:,:,:) :: mapx,mapy


  ! -------------------------------------------------------
  ! blockx by blocky variables
  ! -------------------------------------------------------
  ! Map from SC block to reference block number
  integer, allocatable, dimension(:,:) :: SC_to_refx, SC_to_refy
  ! Map from reference block to SC block number
  integer, allocatable, dimension(:,:) :: ref_to_SCx, ref_to_SCy
  ! Who owns which block (linear processor number), for each process-group
  integer, allocatable, dimension(:,:,:) :: proc_block


  ! -------------------------------------------------------
  ! block_size by block variables
  ! -------------------------------------------------------
  ! Atom numbers in given ref row/column block
  type(block_atom), allocatable, dimension(:,:) :: ref_row_block_atom, &
       ref_col_block_atom
  ! Atom numbers in given SC row/column block
  type(block_atom), allocatable, dimension(:,:) :: SC_row_block_atom, SC_col_block_atom

!!***

contains

    
!!****f* ScalapackFormat/allocate_arrays *
!!
!!  NAME 
!!   allocate_arrays
!!  USAGE
!!   allocate_arrays
!!  PURPOSE
!!   Allocates storage for Scalapack format arrays
!!
!!   Let's be clear about storage and format of matrices in Fortran:
!!   when you go access A(2,3), you go down to the second row, and 
!!   the third entry in that row (i.e. the first index goes down the
!!   column).  I learnt to code mainly in C, so I always have to work
!!   this out - I figured that it might be helpful to have here !
!!
!!   So, when we're allocating stuff below we need the arrays to have
!!   dimension (maxrow, maxcol) so that it has the right dimensions.  But
!!   we'll find reference to "Fortran doing the columns first while C
!!   does the opposite" - all this means is that Fortran goes DOWN the
!!   elements of each column in turn.
!!  INPUTS
!! 
!! 
!!  USES
!!   global_module
!!   generic_comms
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   28th January 2002
!!  MODIFICATION HISTORY
!!   07/02/2002 dave (Happy Birthday Robert)
!!    Added ROBODoc header, removed extraneous write statements
!!   16/04/2002 dave
!!    Added salutation at start of subroutine
!!   2008/05/19 ast
!!    Added timer
!!   2011/02/13 L.Tong
!!    Added k-point parallelisation
!!   2023/08/02 tsuyoshi
!!    Changes for padding of H and S matrices
!!  SOURCE
!!
  subroutine allocate_arrays (nkp)

    use datatypes, ONLY: double
    use global_module, ONLY: iprint_DM, numprocs
    use GenComms, ONLY: cq_abort, myid
    use numbers, ONLY: RD_ERR

    implicit none

    ! passed variables 
    integer, intent(in) :: nkp  ! total number of k-points, used as input to avoid circular dependency with DiagModule

    ! Local variables
    integer :: stat
    logical :: flag_err

    if(iprint_DM>3.AND.myid==0) write(io_lun,fmt='(10x,i5,a)') myid,' Starting Allocate Arrays'

    ! Calculate maximum numbers of blocks in different directions
    if(flag_padH) then
       blocks_r = ceiling((real(matrix_size,double)-RD_ERR)/block_size_r)
       blocks_c = ceiling((real(matrix_size,double)-RD_ERR)/block_size_c)
       if(blocks_r .ne. blocks_c) call cq_abort("ScalapackFormat: blocks_r /= blocks_c")
       matrix_size_padH = blocks_r * block_size_r
    else
       blocks_r = (matrix_size/block_size_r)
       blocks_c = (matrix_size/block_size_c)
       matrix_size_padH = matrix_size
    endif
    if(myid==0.AND.iprint_DM>1) write(io_lun,*) "matrix_size & matrix_size_padH = ",matrix_size, matrix_size_padH
    if(myid==0.AND.iprint_DM>3) write(io_lun,1) blocks_r,blocks_c
    maxrow = floor(real(blocks_r/proc_rows))+1
    maxcol = floor(real(blocks_c/proc_cols))+1
    if(iprint_DM>3.AND.myid==0) write(io_lun,fmt='(10x,a,2i5)') 'maxrow, maxcol: ',maxrow,maxcol
    call start_timer(tmr_std_allocation)
    allocate(mapx(numprocs,maxrow,maxcol),mapy(numprocs,maxrow,maxcol),STAT=stat)
    if(stat/=0) call cq_abort("ScalapackFormat: Could not alloc map",stat)
    allocate(procid(proc_groups, proc_rows,proc_cols),STAT=stat)
    if(stat/=0) call cq_abort("ScalapackFormat: Could not alloc procid",stat)
    allocate(SC_to_refx(blocks_r,blocks_c),SC_to_refy(blocks_r,blocks_c),&
         ref_to_SCx(blocks_r,blocks_c),ref_to_SCy(blocks_r,blocks_c),&
         proc_block(proc_groups,blocks_r,blocks_c),STAT=stat)
    if(stat/=0) call cq_abort("ScalapackFormat: Could not alloc bxb var",stat)
    allocate(ref_row_block_atom(block_size_r,blocks_r),&
         ref_col_block_atom(block_size_c,blocks_c),&
         SC_row_block_atom(block_size_r,blocks_r),&
         SC_col_block_atom(block_size_c,blocks_c),STAT=stat)
    if(stat/=0) call cq_abort("ScalapackFormat: Could not alloc bxa var",stat)
    allocate(CQ2SC_row_info(matrix_size_padH), my_row(matrix_size),proc_start(numprocs), STAT=stat)
    nprocs_max = floor(real(numprocs/proc_groups))+1
    nkpoints_max = floor(real(nkp/proc_groups))+1
    allocate(pg_procs(proc_groups,nprocs_max),pg_kpoints(proc_groups,nkpoints_max),STAT=stat)
    if(stat/=0) call cq_abort("ScalapackFormat: Could not alloc pg_procs, pg_kpoints",stat)
    allocate(N_procs_in_pg(proc_groups),N_kpoints_in_pg(proc_groups),STAT=stat)
    if(stat/=0) call cq_abort("ScalapackFormat: Could not alloc N_procs_in_pg, N_kpoints_pg",stat)
    allocate(pgroup(numprocs),STAT=stat)
    if(stat/=0) call cq_abort("ScalapackFormat: Could not alloc pgroup",stat)
    call stop_timer(tmr_std_allocation)
    return
1   format(10x,'AllocArr: block sizes are: ',2i5)
  end subroutine allocate_arrays
!!***

!!****f* ScalapackFormat/deallocate_arrays *
!!
!!  NAME 
!!   deallocate_arrays
!!  USAGE
!!   deallocate_arrays
!!  PURPOSE
!!   Deallocates the memory allocated by allocate_arrays
!!  INPUTS
!! 
!!  USES
!!   global_module, GenComms
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   April/May 2002
!!  MODIFICATION HISTORY
!!   17/06/2002 dave
!!    Added ROBODoc header
!!   2008/05/19 ast
!!    Added timer
!!   2011/02/13 L.Tong
!!    Added k-point parallelisation
!!  SOURCE
!!
  subroutine deallocate_arrays

    use global_module, ONLY: iprint_DM
    use GenComms, ONLY: cq_abort, myid

    implicit none

    ! Local variables
    integer :: stat

    if(iprint_DM>3.AND.myid==0) write(io_lun,fmt='(10x,a)') myid,' Starting Deallocate Arrays'
    call start_timer(tmr_std_allocation)
    deallocate(CC_to_SC,CQ2SC_row_info, my_row,proc_start, STAT=stat)
    if(stat/=0) call cq_abort("ScalapackFormat: Could not dealloc CC2SC, CQ2SC",stat)
    deallocate(ref_row_block_atom,&
         ref_col_block_atom,&
         SC_row_block_atom,&
         SC_col_block_atom,STAT=stat)
    if(stat/=0) call cq_abort("ScalapackFormat: Could not dealloc bxa var",stat)
    deallocate(SC_to_refx,SC_to_refy,&
         ref_to_SCx,ref_to_SCy,&
         proc_block,STAT=stat)
    if(stat/=0) call cq_abort("ScalapackFormat: Could not dealloc bxb var",stat)
    deallocate(procid,STAT=stat)
    if(stat/=0) call cq_abort("ScalapackFormat: Could not dealloc procid",stat)
    deallocate(mapx,mapy,STAT=stat)
    if(stat/=0) call cq_abort("ScalapackFormat: Could not dealloc map",stat)
    deallocate(N_procs_in_pg,N_kpoints_in_pg,STAT=stat)
    if(stat/=0) call cq_abort("ScalapackFormat: Could not dealloc N_procs_in_pg, N_kpoints_in_pg",stat)
    deallocate(pg_procs,pg_kpoints,STAT=stat)
    if(stat/=0) call cq_abort("ScalapackFormat: Could not dealloc pg_procs, pg_kpoints",stat)
    deallocate(pgroup,STAT=stat)
    if(stat/=0) call cq_abort("ScalapackFormat: Could not dealloc pgroup",stat)
    call stop_timer(tmr_std_allocation)
    return
  end subroutine deallocate_arrays
!!***

!!****f* ScalapackFormat/pg_initialise *
!!
!!  NAME
!!   pg_initialise
!!  USAGE 
!!   pg_initialise
!!  PURPOSE
!!   Distribute the processes and kpoints into proc_groups. 
!!  INPUTS
!!
!!  USES
!!  
!!  AUTHOUR
!!   L.Tong
!!  CREATION DATE 
!!   13th Feburary 2011
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine pg_initialise (nkp)
    
    use GenComms, ONLY: cq_abort, inode
    use global_module, ONLY: numprocs
  
    implicit none
    
    ! passed variables
    integer, intent(in) :: nkp  ! number of k-points, use as input to avoid circular dependencies with DiagModule
    
    ! local variables
    integer :: stat
    integer :: ng, np, nk
    integer, dimension(:), allocatable :: counter
    
    pgid = mod ((inode - 1), proc_groups) + 1  ! get the local proc_group id
  
    ! temporary limitation, the process-groups must have same number of processes
    if (mod (numprocs, proc_groups) > 0) &
         call cq_abort ('pg_initialise: the number of processes in each group must be the same.', mod (numprocs, proc_groups)) 

    do ng = 1, proc_groups
       ! calculate the total number of members in proc_group ng
       if (ng <= mod (numprocs, proc_groups)) then
          N_procs_in_pg(ng) = nprocs_max  ! nprocs_max is worked out in subroutine allocate_array
       else
          N_procs_in_pg(ng) = nprocs_max - 1
       end if
       ! calculate the total number of kpoints assigned to proc_group ng
       if (ng <= mod (nkp, proc_groups)) then
          N_kpoints_in_pg(ng) = nkpoints_max
       else
          N_kpoints_in_pg(ng) = nkpoints_max - 1
       end if
    end do
    
    ! temporary limitation: the total number of processes in each group must equal to proc_rows*proc_cols
    if (N_procs_in_pg(pgid) /= proc_rows*proc_cols) &
         call cq_abort ('pg_initialise: Number of nodes in each proc_group must equal to proc_rows*proc_cols', &
         N_procs_in_pg(pgid), proc_rows*proc_cols)
    
    allocate (counter(proc_groups), STAT=stat)
    if (stat /= 0) call cq_abort ("ScalapackFormat, deallocate_arryays: Failed to allocate counter",stat)
    counter = 1
    pg_procs = 0
    do np = 1, numprocs 
       pgroup(np) = mod ((np - 1), proc_groups) + 1
       ng = pgroup(np)
       pg_procs(ng, counter(ng)) = np
       counter(ng) = counter(ng) + 1
    end do

    counter = 1
    pg_kpoints = 0
    do nk = 1, nkp
       ng = mod ((nk - 1), proc_groups) + 1
       pg_kpoints(ng, counter(ng)) = nk
       counter(ng) = counter(ng) + 1
    end do
    deallocate (counter, STAT=stat)
    if (stat /= 0) call cq_abort ("ScalapackFormat, deallocate_arryays: Failed to deallocate counter",stat)
    return
  end subroutine pg_initialise

!!****f* ScalapackFormat/ref_to_SC_blocks *
!!
!!  NAME 
!!   ref_to_SC_blocks
!!  USAGE
!!   ref_to_SC_blocks
!!  PURPOSE
!!   Takes blocks allocated to processor grid in reference format 
!!   and converts to compacted Scalapack blocks
!!  INPUTS
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   28th January 2002
!!  MODIFICATION HISTORY
!!   07/02/2002 dave (Happy Birthday Robert)
!!    Generalised to different row/col blocks sizes
!!   14/03/2002 dave
!!    Added better output for different iprint levels
!!   16/04/2002 dave
!!    Added salutation at start of subroutine
!!   01/05/2002 dave
!!    Changed so that Scalapack grid is only written by one processor
!!   2008/05/19 ast
!!    Added timer
!!  SOURCE
!!
  subroutine ref_to_SC_blocks 

    use GenComms, ONLY: cq_abort, myid
    use global_module, ONLY: iprint_DM

    implicit none

    integer :: i, j, n, nrow, ncol, prow, pcol, proc
    integer :: ng

    if(iprint_DM>3.AND.myid==0) write(io_lun,fmt='(8x,a)') 'Starting Ref To SC Blocks'
    ! Construct processor ids
    if(iprint_DM>2.AND.myid==0) write(io_lun,fmt="(8x,'Scalapack Processor Grid')") 
    do ng = 1, proc_groups
       n = 1
       do i=1,proc_rows
          do j=1,proc_cols
             procid(ng, i, j) = pg_procs(ng, n)
             if(n>N_procs_in_pg(ng)) call cq_abort('Ref2SC: Too many processors in group',ng,n)
             n = n + 1
          end do
          if(iprint_DM>4.AND.myid==0) write(io_lun,*) (procid(ng,i,j),j=1,proc_cols)
       end do
    end do
    ! now build list of blocks and where they go
    if(iprint_DM>3.AND.myid==0) &
         write(io_lun,fmt="(8x,'Map from local block on processor (first three) to reference block (last two)')")
    if(iprint_DM>3.AND.myid==0) write(io_lun,fmt="(27x,'Proc Nrow Ncol Mapx Mapy')") 
    do ng = 1, proc_groups
       do i=1,blocks_r                          ! Rows of blocks in ref format
          prow = mod(i-1,proc_rows)+1             ! Processor row for this row
          nrow = floor(real((i-1)/proc_rows))+1    ! Which row on the processor 
          do j=1,blocks_c                       ! Cols of blocks in ref format
             pcol = mod(j-1,proc_cols)+1          ! Processor col for this row
             ncol = floor(real((j-1)/proc_cols))+1 ! Which col on the processor
             ! Which processor (linear number) is this ?
             proc = procid(ng, prow, pcol)
             ! Map from block on processor to reference format
             mapx(proc,nrow,ncol) = i
             mapy(proc,nrow,ncol) = j
             proc_block(ng, i, j) = proc ! Owner of block (linear number)
             if(iprint_DM>3.AND.myid==0) &
                  write(io_lun,fmt="(8x,'Proc: ',i5,'Mapx and y: ',6i5)") myid,ng,proc,nrow,ncol,i,j
          end do
       end do
    end do
    return
  end subroutine ref_to_SC_blocks
!!***

!!****f* ScalapackFormat/make_maps *
!!
!!  NAME 
!!   make_maps
!!  USAGE
!!   make_maps
!!  PURPOSE
!!   Creates maps between blocks on processors and reference 
!!   format.  Loops over processor rows, and blocks in that
!!   processor's chunk, then over processor cols and blocks in
!!   that chunk (again), and uses the fact that as you do this
!!   you'll be traversing the Scalapack columns in order (and
!!   then rows).  This diagram may help:
!!      c1 2 3 4 5
!!
!!   r1  1 1 2 2 3
!!   r2  1 1 2 2 3
!!   r3  1 1 2 2 3
!!   r4  4 4 5 5 6
!!   r5  4 4 5 5 6
!! 
!!   This is a matrix with 5 blocks in each row and column, 
!!   and a (2x3) processor grid. With processor row set to 1, 
!!   and row set to 1, then we loop over processor columns from
!!   1 to 3, and cols on each processor, we get the following 
!!   pattern, for row r1:
!! 
!!   proc col = 1, cols = 1->2 gives us c1, c2 with m=1,2
!!              2, cols = 1->2 gives us c3, c4 with m=3,4
!!              3, cols = 1    gives us c5     with m=5
!!   This repeats for rows 2 and 3 on the processor (r2, r3).  Then 
!!   we increment processor row, and for rows 1 and 2 get r4 and r5, 
!!   with the columns repeating as above.
!!  INPUTS
!!   None
!!  USES
!!   global_module
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   28th January 2002
!!  MODIFICATION HISTORY
!!   07/02/2002 dave (Happy Birthday Robert)
!!    Generalised to different row/col blocks sizes
!!   07/03/2002 dave
!!    Changed so that we go over row procs and blocks and then
!!    col procs and blocks.  I think that this makes more sense.
!!    Added my_row assignments (17.10)
!!   15/03/2002 dave
!!    Significantly increased description of purpose, removed
!!    a few unnecessary comments, tidied and rearranged and
!!    added comments.
!!   2008/05/19 ast
!!    Added timer
!!   2011/02/13 L.Tong
!!    Added modification required for k-point parallelisation
!!   2013/04/24 15:07 dave
!!    Bug fix (initially due to Conn O'Rourke) to degine proc properly
!!  SOURCE
!!
  subroutine make_maps 

    use global_module, ONLY: iprint_DM
    use GenComms, ONLY: myid

    implicit none

    integer :: i,j,m,n,row,col,proc,ng
    integer :: row_max_n, col_max_n, loc_max_row, loc_max_col

    if(iprint_DM>3.AND.myid==0) write(io_lun,3)
    if(iprint_DM>3.AND.myid==0) write(io_lun,4)
    ! The first row_max_n procs have 1 more block than the rest
    row_max_n = mod(blocks_r,proc_rows)
    col_max_n = mod(blocks_c,proc_cols)
    if(iprint_DM>3.AND.myid==0) write(io_lun,fmt='(8x,a,2i4)') 'N for row, col: ',row_max_n, col_max_n
    if(iprint_DM>3.AND.myid==0) write(io_lun,fmt='(8x,a,2f4.0)') 'Loc_max_row, col: ',&
         aint(real(blocks_r/proc_rows)),aint(real(blocks_c/proc_cols))
    my_row = 0
    ! first record proc_start(:)%rows and %cols map, this is proc proc dependent
    ! n is block row in global SC format
    n = 1
    do i = 1, proc_rows ! looping over proc-grid rows
       if(i<=row_max_n) then
          loc_max_row = floor(real(blocks_r/proc_rows))+1
       else
          loc_max_row = floor(real(blocks_r/proc_rows))
       end if
       do row = 1, loc_max_row ! looping over local SC format block rows
          do j = 1, proc_cols ! looping over proc-grid cols
             if (j <= col_max_n) then
                loc_max_col = floor(real(blocks_c/proc_cols))+1
             else
                loc_max_col = floor(real(blocks_c/proc_cols))
             end if
             do ng = 1, proc_groups
                proc = procid(ng, i, j)
                if(proc == myid + 1) my_row(n) = 1  ! This is for receiving data
                ! Record rows and columns on this processor
                proc_start(proc)%rows = loc_max_row
                proc_start(proc)%cols = loc_max_col
             end do
          end do
          if (iprint_DM>3.AND.myid==0) write(io_lun,2) myid,n,blocks_r,my_row(n)
          n = n + 1 ! accumulate global SC block row counter
          if (n > blocks_r) n = 1       
       end do
    end do
    
    ! now make the SC_to_ref and ref_to_SC maps that are identical to all groups, 
    ! hence only needs to work on one group
    ! n and m are row and column block in SC format, reset
    n = 1
    m = 1
    ng = 1  ! using the first proc_group
    do i = 1, proc_rows ! looping over proc-grid
       proc = procid(ng, i, 1) ! get global proc id
       do row = 1, proc_start(proc)%rows  ! looping over local SC block row *The number of rows is independent of proc_cols*
          do j = 1, proc_cols
             proc = procid(ng, i, j) ! get global proc id
             do col = 1, proc_start(proc)%cols  ! looping over local SC block col
                ! For this SC block, record reference block coordinates
                SC_to_refx(n,m) = mapx(proc,row,col)
                SC_to_refy(n,m) = mapy(proc,row,col)
                ! Use the mapping to record SC block coords for ref block
                ref_to_SCx(mapx(proc,row,col),mapy(proc,row,col)) = n
                ref_to_SCy(mapx(proc,row,col),mapy(proc,row,col)) = m
                ! Write out if necessary
                if(iprint_DM>3.AND.myid==0) write(io_lun,1) proc,row,col,mapx(proc,row,col),mapy(proc,row,col),n,m
                m = m + 1
                if(m>blocks_c) m=1 ! this is need so that the counter restarts for every block row
             end do
          end do
          n = n + 1
          if(n>blocks_r) n=1             
       end do
    end do
    return
1   format(8x,'MakeMaps: proc,row,col,map,n,m: ',7i5)
2   format(8x,'MakeMaps: proc,n,max,my_row: ',4i5)
3   format(8x,'Welcome to make_maps')
4   format(8x,'Here we will create maps from ref to SC and back again')
  end subroutine make_maps
!!***

!!****f* ScalapackFormat/find_SC_row_atoms *
!!
!!  NAME 
!!   find_SC_row_atoms
!!  USAGE
!!   find_SC_row_atoms
!!  PURPOSE
!!   Distributes atoms over rows in Scalapack format according
!!   to Conquest processor/primary atom ordering (i.e. loop over
!!   processors, and then over primary atoms for each processor,
!!   and that's the order of the rows in Scalapack format).  We
!!   insist on this ordering to make data distribution from CQ
!!   compressed row format to Scalapack format easier - it'll 
!!   involve reordering of columns on-node, but no more.
!!
!!   Let's think about how to get Conquest ordering of atoms. The global
!!   atom numbering is given by id_glob (but the index is CC ?). Within
!!   primary_module, we go over atoms in the primary set partition by
!!   partition (using groups%ng_on_node), and then find the CC label of
!!   the partition using groups%ngnode.  We then use id_glob and x,y,z to
!!   identify the atom.  What will we need to identify the atom ? Processor, 
!!   partition, atom in partition would all be useful.  Actually, if we have
!!   the CC label of the partition then we have the processor and sequence 
!!   number of the partition on the processor.
!!
!!   Looping over processors and partitions on processor, and atoms in the
!!   partition will give us the correct order.  The structure groups is going
!!   to be vital here - groups%ng_on_node tells us how many partitions there
!!   are, while groups%ngnode gives the CC label.  Then groups%nm_group is
!!   the number of atoms, and this will be the desired order.
!!
!!   The next question is what we store - I think that the CC label of the
!!   partition (which gives us processor and sequence no) and sequence number
!!   of atom in partition should be OK.
!!
!!   Note that we loop over support functions on each atom (which, of course,
!!   correspond to rows of the matrix) when we consider the distribution of
!!   data over Scalapack blocks.  This means that data from any given atom may
!!   be spread across different blocks, and that atoms are free to have varying
!!   numbers of support functions.
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   28th January 2002
!!  MODIFICATION HISTORY
!!   07/02/2002 dave (Happy Birthday Robert)
!!    Generalised to different row/col blocks sizes
!!   19/02/2002 dave (Happy Birthday Christopher)
!!    Added Conquest/ParaDens atom numbering
!!   20/02/2002 dave
!!    Added structure members to SC_row_block_atom
!!   19/03/2002 dave
!!    Rewrote entirely - first, build CQ2SC_row_info which gives all information
!!    we could want about the Conquest atom corresponding to this row.  Then from
!!    this build up SC_row_block_atom, working by rows (not atoms as before). This
!!    frees us to have matrix size equal to \sum_i nsf(i), where nsf(i) is the 
!!    number of support functions on atom i
!!   16/04/2002 dave
!!    Added salutation at start of subroutine
!!   2008/05/19 ast
!!    Added timer
!!   2023/07/20 tsuyoshi
!!    Change for padding H and S matrices  (to improve the efficiency of diagonalization)
!!  SOURCE
!!
  subroutine find_SC_row_atoms 

    use gencomms, ONLY: cq_abort, myid
    use group_module, ONLY: parts
    use global_module, ONLY: numprocs, iprint_DM
    use species_module, ONLY: nsf_species, species

    implicit none

    ! Local variables
    integer :: i, patom, part, proc, CC, brow, SCblock, supfn
    integer :: num_elem_pad, nproc_pad, iblock_pad, end_elem_pad, start_elem_pad

    ! For padding H and S matrices
    num_elem_pad = matrix_size_padH - matrix_size     ! # of padded elements
    nproc_pad = mod( blocks_r-1, proc_rows ) + 1      ! ID of the responsible process for the padded part
    iblock_pad = sum( proc_start( 1 : nproc_pad )%rows ) ! # of blocks before reaching the block including the padded part
    end_elem_pad = block_size_r * iblock_pad          ! element ID for the end of the padded part
    start_elem_pad = end_elem_pad - num_elem_pad + 1  ! element ID for the start of the padded part

    if(iprint_DM>3.AND.myid==0) write(io_lun,fmt='(10x,i5,a)') myid,' Starting Find SC Row Atoms'
    ! -----------------------------------------------------------------
    ! Loop over matrix using processor/partition/sequence order
    i = 1 ! Indexes matrix row
    ! For padding H and S matrices:  
    !  CQ2SC_row_info(:)%CClabel is now initialised and some of the values can be kept 0 for the padded parts.
    !  Other members are basically not used for padded parts, but also initialised to be 0 for safety.
     CQ2SC_row_info(:)%CClabel = 0
     CQ2SC_row_info(:)%support_fn = 0
     CQ2SC_row_info(:)%atom = 0
     CQ2SC_row_info(:)%partition = 0
    
    do proc=1,numprocs
       proc_start(proc)%startrow = i ! This is where data for the processor starts
       if(parts%ng_on_node(proc)>0) then
          do part=1,parts%ng_on_node(proc)  ! Loop over partitions
             CC = parts%ngnode(parts%inode_beg(proc)+part-1)
             if(parts%nm_group(CC)>0)then
                do patom=1,parts%nm_group(CC) ! Loop over atoms in partition
                   do supfn=1,nsf_species(species(parts%icell_beg(CC)+patom-1))!nsf ! Generalise to atom-dependent sup fn number !
                      CQ2SC_row_info(i)%support_fn = supfn
                      CQ2SC_row_info(i)%atom = patom
                      CQ2SC_row_info(i)%partition = part
                      CQ2SC_row_info(i)%CClabel = CC
                      i = i+1
                      if( i==start_elem_pad ) i=i+num_elem_pad
                      if(i>matrix_size_padH+1) call cq_abort('Too many support functions !',i)
                   end do
                end do
             end if
          end do
       end if
    end do
    ! End loop over matrix using proc/part/seq
    ! -----------------------------------------------------------------
    ! Now loop over matrix using blocks (same order, different indices)
    i=1 ! Indexes matrix row
    proc = 1
    do SCblock=1,blocks_r ! Blocks
       do brow=1,block_size_r ! Rows in block
          if(i==proc_start(proc)%startrow) then  ! Note when a processor's data starts
             proc_start(proc)%block = SCblock
             proc_start(proc)%atom  = brow
             proc=proc+1
             if(proc>numprocs) proc = numprocs
          end if
          ! Record SupFn/Atom/Partition for this row and block
          SC_row_block_atom(brow,SCblock)%part  = CQ2SC_row_info(i)%CClabel
          SC_row_block_atom(brow,SCblock)%atom  = CQ2SC_row_info(i)%atom
          SC_row_block_atom(brow,SCblock)%support_fn = CQ2SC_row_info(i)%support_fn
          if(iprint_DM>3.AND.myid==0) &
               write(io_lun,2) myid,CQ2SC_row_info(i)%CClabel,CQ2SC_row_info(i)%atom,CQ2SC_row_info(i)%support_fn
          i=i+1
          if(i>matrix_size_padH+1) call cq_abort('Too many support functions !',i)
       end do
    end do
    ! End loop over matrix using blocks
    ! -----------------------------------------------------------------
    return

2   format(8x,'Proc: ',i5,' Part, Seq, SupFn: ',3i5)
  end subroutine find_SC_row_atoms
!!***

!!****f* ScalapackFormat/find_ref_row_atoms *
!!
!!  NAME 
!!   find_ref_row_atoms
!!  USAGE
!!   find_ref_row_atoms 
!!  PURPOSE
!!   Given the Scalapack format listing for atoms corresponding
!!   to each row, generate the reference format listing for the
!!   atoms corresponding to each row, and then use the symmetry
!!   of the reference format to find the atoms corresponding to
!!   each column in reference format.
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   28th January 2002
!!  MODIFICATION HISTORY
!!   07/02/2002 dave (Happy Birthday Robert)
!!    Generalised to different row/col blocks sizes
!!   08/02/2002 dave
!!    Fixed assumption about block size symmetry
!!   08/02/2002 dave
!!    Changed names of conversion arrays and added comments
!!   20/02/2002 dave
!!    Added structure members to ref_row_block_atom
!!   04/03/2002 dave
!!    Changed CC_to_SC map to allow row storage as well as column
!!   05/04/2002 dave
!!    Many changes over last month: mainly changed so that we work
!!    row-by-row (i.e. by support function) rather than atom-by-atom
!!   16/04/2002 dave
!!    Added salutation at start of subroutine
!!   2008/05/19 ast
!!    Added timer
!!   2013/07/05 dave
!!    Moved check on cb overrunning so that erroneous error goes away
!!   2023/7/20 tsuyoshi
!!    Changed for the version of padding H matrix
!!  SOURCE
!!
  subroutine find_ref_row_atoms 

    use global_module, ONLY: iprint_DM
    use gencomms, ONLY: cq_abort, myid

    implicit none

    ! Local variables
    integer :: rb, cb, SCblockx, blockrow, blockcol
    integer :: part, seq, supfn

    if(iprint_DM>3.AND.myid==0) write(io_lun,fmt='(10x,i5,a)') myid,' Starting Find Ref Row Atoms'
    blockcol = 1
    cb = 1
    ! For padding of H matrix, we need to initialize ref_row_block_atom and ref_col_block_atom,
    ! because they need to be kept 0 for the padding part.  
    !  (it might be safer to initialize other members (atom and support_fn)
    ref_row_block_atom(:,:)%part = 0
    ref_col_block_atom(:,:)%part = 0
    ! for safety, other members are also initialised.
    ref_row_block_atom(:,:)%atom = 0
    ref_col_block_atom(:,:)%atom = 0
    ref_row_block_atom(:,:)%support_fn = 0
    ref_col_block_atom(:,:)%support_fn = 0

    if(iprint_DM>3.AND.myid==0) write(io_lun,fmt='(10x,a,2i5)') '  blocks, size: ',blocks_r, block_size_r
    ! Loop over reference row blocks
    do rb = 1,blocks_r
       SCblockx = ref_to_SCx(rb,1)  ! find equivalent block in SC format
       if(iprint_DM>3.AND.myid==0) write(io_lun,2) myid,rb,SCblockx
       do blockrow = 1,block_size_r ! Loop over rows in block
          part = SC_row_block_atom(blockrow,SCblockx)%part
          !For the padded part, setting of ref_row_block_atom should not be done. (kept to be 0)
          if(part == 0) cycle   
          seq  = SC_row_block_atom(blockrow,SCblockx)%atom
          supfn  = SC_row_block_atom(blockrow,SCblockx)%support_fn
          if(cb>blocks_c+1.OR.cb==blocks_c+1.AND.blockcol>1) then
             call cq_abort("ScalapackFormat: Too many column blocks")
          end if
          ref_row_block_atom(blockrow,rb)%part = part
          ref_row_block_atom(blockrow,rb)%atom = seq
          ref_row_block_atom(blockrow,rb)%support_fn = supfn
          ref_col_block_atom(blockcol,cb)%part = part
          ref_col_block_atom(blockcol,cb)%atom = seq
          ref_col_block_atom(blockcol,cb)%support_fn = supfn
          CC_to_SC(part,seq,supfn)%block_r = SCblockx
          CC_to_SC(part,seq,supfn)%row_r = blockrow
          if(iprint_DM>3.AND.myid==0) write(io_lun,1) myid,part,seq,supfn,rb,blockrow
          blockcol = blockcol + 1
          if(blockcol>block_size_c) then ! End of block - increment
             cb = cb + 1
             blockcol = 1
          end if
       end do
    end do
1   format(8x,'RefRowAt Proc: ',i5,' Part, seq, SupFn: ',3i5,' Block, row: ',2i5)
2   format(8x,'RefRowAt Proc: ',i5,' Ref block, SC block: ',2i5)
    return
  end subroutine find_ref_row_atoms
!!***

!!****f* ScalapackFormat/find_SC_col_atoms *
!!
!!  NAME 
!!   find_SC_col_atoms  
!!  USAGE
!!   find_SC_col_atoms  
!!  PURPOSE
!!   Given the atoms corresponding to the columns in the 
!!   reference format, generate the atoms corresponding to
!!   the columns in Scalapack format.
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   28th January 2002
!!  MODIFICATION HISTORY
!!   07/02/2002 dave (Happy Birthday Robert)
!!    Generalised to different row/col blocks sizes
!!   08/02/2002 dave
!!    Changed names of conversion arrays and added comments
!!   20/02/2002 dave
!!    Added structure members to SC_col_block_atom
!!   28/02/2002 dave
!!    Created mapping from Conquest CC labels to SC col number
!!   04/03/2002 dave
!!    Changed CC_to_SC map
!!   05/04/2002 dave
!!    Many changes over last month: mainly related to conversion from
!!    atom-by-atom to row-by-row (i.e. by individual support functions)
!!   16/04/2002 dave
!!    Added salutation at start of subroutine
!!   2008/05/19 ast
!!    Added timer
!!   2023/07/20 tsuyoshi
!!    Change for padding Hamiltonian and overlap matrices
!!  SOURCE
!!
  subroutine find_SC_col_atoms  

    use global_module, ONLY: iprint_DM
    use group_module, ONLY: parts
    use GenComms, ONLY: myid
    use species_module, ONLY: nsf_species, species
    
    implicit none

    ! Local variables
    integer :: cb, blockcol, refc, part, seq, supfn, np_in_cell

    ! for padding H matrix, we need to initialise SC_col_block_atom
          SC_col_block_atom(:,:)%part = 0
     ! for safety other members are also initialised.
          SC_col_block_atom(:,:)%atom = 0
          SC_col_block_atom(:,:)%support_fn = 0
         
    if(iprint_DM>3.AND.myid==0) write(io_lun,fmt='(10x,i5, a)') myid,' Starting Find SC Col Atoms'
    ! Loop over SC blocks
    do cb = 1,blocks_c
       refc = SC_to_refy(1,cb) ! find equivalent number in reference format
       if(iprint_DM>3.AND.myid==0) write(io_lun,2) myid, cb, refc
       do blockcol = 1,block_size_c
          part = ref_col_block_atom(blockcol,refc)%part  
          ! for the padded part (padding Hmatrix version)
          if(part == 0) cycle    ! for padding H and S matrices  (for padded part)
          seq  = ref_col_block_atom(blockcol,refc)%atom  
          supfn  = ref_col_block_atom(blockcol,refc)%support_fn
          SC_col_block_atom(blockcol,cb)%part = part
          SC_col_block_atom(blockcol,cb)%atom = seq
          SC_col_block_atom(blockcol,cb)%support_fn = supfn
          CC_to_SC(part,seq,supfn)%block_c = cb
          CC_to_SC(part,seq,supfn)%row_c  = blockcol
       end do
    end do
    if(iprint_DM>3.AND.myid==0) then
       np_in_cell = parts%ngcellx*parts%ngcelly*parts%ngcellz
       write(io_lun,4) myid
       do part=1,np_in_cell  ! Loop over partitions
          do seq=1,parts%nm_group(part) ! Loop over atoms in partition
             do supfn=1,nsf_species(species(parts%icell_beg(part)+seq-1))
                write(io_lun,3) part,seq,supfn,CC_to_SC(part,seq,supfn)
             end do
          end do
       end do
    endif
    return

2   format(8x,'Proc: ',i5,' SC col, ref col: ',2i5)
3   format(8x,3i5,' : ',4i5)
4   format(8x,'On proc ',i3,' CC_to_SC is',/,2x,' Part Atom  SFn   Label')
  end subroutine find_SC_col_atoms
!!***
end module ScalapackFormat
