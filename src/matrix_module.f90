! $Id$
! -----------------------------------------------------------
! Module matrix_module
! -----------------------------------------------------------
! Code area 2: matrices
! -----------------------------------------------------------

!!****h* Conquest/matrix_module
!!  NAME
!!   matrix_module
!!  PURPOSE
!!   This module defines various derived types associated with matrices (matrix
!!   for generic indexing, matrix_trans for variables associated with the local
!!   transpose and matrix_halo for variables associated with the halo of a 
!!   matrix).  It also includes subroutines associated with these types for 
!!   allocating and deallocating memory
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   24/11/99
!!  MODIFICATION HISTORY
!!   03/02/00 by D.R.Bowler
!!    Now all data in type(matrix) is stored by partition for efficient data
!!    transfer
!!   12/04/00 by DRB
!!    Added itran_addr to trans_remote type
!!   21/06/2001 dave
!!    Added ROBODoc headers throughout and changed to use cq_abort
!!   14:38, 26/02/2003 drb 
!!    Tidied cq_aborts and removed integ in halo structure
!!   07:59, 2003/06/11 dave
!!    Added TM's debugging changes
!!   2008/02/06 08:23 dave
!!    Changed for output to file not stdout
!!   2008/05/22 ast
!!    Added timers
!!   2012/01/18 16:52 dave
!!    Added blip transfer derived type
!!   2014/09/15 18:30 lat
!!    fixed call start/stop_timer to timer_module (not timer_stdlocks_module !)
!!***
module matrix_module

  ! Module usage
  use datatypes
  use basic_types
  use global_module,          only: io_lun
  use GenComms,               only: cq_abort
  use timer_module,           only: start_timer,stop_timer
  use timer_stdclocks_module, only: tmr_std_allocation

  implicit none

!!****s* matrix_module/matrix *
!!  NAME
!!   matrix
!!  PURPOSE
!!   Combines together all variables needed to index matrices
!!  AUTHOR
!!   D.R.Bowler
!!  SOURCE
!!
  type matrix
     integer :: mx_nab,mx_abs
     integer :: sf1_type, sf2_type
     integer :: length
     real(double),pointer :: radius(:)
     integer(integ) :: array_posn ! Starting position in array
     integer(integ) :: offset ! This plus i_acc gives us i_dir
     integer(integ) :: nd_offset ! This plus i_acc gives us i_dir
     integer(integ) :: n_atoms, part_nabs, part_nd_nabs ! Lengths of n_nab, i_acc, i_nd_acc
     integer(integ),pointer :: i_part(:) ! Allocated
     integer(integ),pointer :: onsite(:)
     ! From here on, we point to the comms array
     integer(integ),pointer :: n_nab(:)  ! Neighbours of primary atom i
     integer(integ),pointer :: ndimi(:)  ! nsf for neighbours
     integer(integ),pointer :: ndimj(:)  ! nsf for neighbours
     integer(integ),pointer :: i_seq(:)  ! part-seq no. of atom in neigh list
     integer(integ),pointer :: i_acc(:)  ! accumulator for atom in neigh list
     integer(integ),pointer :: i_nd_acc(:)  ! accumulator for SFs in neigh list
     integer(integ),pointer :: npxyz(:)  ! offset of CS part of atom in n. list
  end type matrix
!!***

!!****s* matrix_module/matrix_trans *
!!  NAME
!!   matrix_trans
!!  PURPOSE
!!   Combines variables needed to perform and index local transpose
!!
!!   N.B. iprim is indexed as iprim(nn), nn=i_beg(ia)+i-1, with ia the
!!   halo sequence number and i the neighbour of ia
!!  AUTHOR
!!   D.R.Bowler
!!  SOURCE
!!
  type matrix_trans
     integer(integ),pointer :: i_prim(:)  ! Prim atom spec as neig of halo atom
     integer(integ),pointer :: i_beg(:)   ! Accum of n_hnab
     integer(integ),pointer :: i_nd_beg(:)   ! Accum of n_hnab
     integer(integ),pointer :: n_hnab(:)  ! No of prim atoms neigh halo atom
     ! This type's values for maxima
     integer :: mx_halo,mx_nab
  end type matrix_trans
!!***

!!****s* matrix_module/matrix_halo *
!!  NAME
!!   matrix_halo
!!  PURPOSE
!!   Collects all variables and indices associated with a halo
!!
!!   A halo is that collection of atoms which are within range of 
!!   at least one atom in the primary set.  Halo partitions contain 
!!   at least one halo atom, and halo nodes are responsible for at 
!!   least one halo partition
!!  AUTHOR
!!   D.R.Bowler
!!  SOURCE
!!
  type matrix_halo
     integer :: np_in_halo ! Partns in halo
     integer :: ni_in_halo ! Atoms in halo
     integer,allocatable :: nh_part(:)    ! No of halo atoms in halo part
     integer,allocatable :: j_beg(:)      ! accumulator of nh_part
     integer,allocatable :: lab_hcell(:)  ! sim cell part (CC) = halo part
     integer,allocatable :: lab_hcover(:) ! CS part (CC) = halo part
     integer,allocatable :: j_seq(:)      ! Part seq of halo atom
     integer,allocatable :: i_h2d(:)      ! halo->neigh trans for halo atom
     integer,allocatable :: i_halo(:)     ! halo seq of atom in CS
     integer,allocatable :: i_hbeg(:)     ! Where CS part starts in i_halo
     integer,allocatable :: ndimi(:)
     integer,allocatable :: ndimj(:)
     ! This type's values for maxima
     integer :: mx_part,mx_halo
  end type matrix_halo
!!***

!!****s* matrix_module/comms_data *
!!  NAME
!!   comms_data
!!  PURPOSE
!!   Communication data for a multiplication
!!  AUTHOR
!!   D.R.Bowler
!!  SOURCE
!!
  type comms_data
     integer, pointer :: np_send(:)         ! No part to send to neig proc
     integer, pointer :: pl_send(:,:)       ! List of these parts
     integer, pointer :: neigh_node_list(:) ! Neig node for part
     integer, pointer :: ilen2rec(:,:)      ! Len of data to recv from proc
     integer, pointer :: ilen3rec(:,:)      ! Len of data to recv from proc
     integer, pointer :: istart(:,:)        ! Start of data to recv
     integer, pointer :: ncomm(:)           ! Real proc no of neig proc
     integer :: inode                       ! Number of neig proc
     integer :: int_win, ele_win  ! These are the windows for MPI2
  end type comms_data
!!***

!!****s* matrix_module/blip_transfer *
!!  NAME
!!   blip_transfer
!!  PURPOSE
!!   Communication data for a blip transfers (used in analytic blip integration)
!!  AUTHOR
!!   D.R.Bowler
!!  SOURCE
!!
  type blip_transfer
     integer :: nproc                   ! Number of processes to communicate with
     integer :: npart_send              ! Total number of partitions to send
     integer, pointer :: ncomm(:)       ! Processes to communicate with
     integer, pointer :: np_send(:)     ! Number of partitions to send
     integer, pointer :: neigh_pl(:)    ! List of processors for partitions
     integer, pointer :: pl_send(:,:)   ! List of partitions to send
     integer, pointer :: partst(:)      ! Start of data for partition in coefficient_array 
     integer, pointer :: partlen(:)     ! Size of data in partition 
     integer, pointer :: len_send(:,:)  ! Length to send
     integer, pointer :: len_recv(:)  ! Length to receive
  end type blip_transfer
!!***

!!****s* matrix_module/matrix_mult *
!!  NAME
!!   matrix_mult
!!  PURPOSE
!!   Holds useful indices for matrix mult
!!  AUTHOR
!!   D.R.Bowler
!!  SOURCE
!!
  type matrix_mult
     integer :: mult_type          ! Extension or reduction (1 or 2)
     type(group_set), pointer :: parts  ! Parts type
     type(primary_set), pointer :: prim ! Primary set type
     type(cover_set), pointer :: gcs    ! CS type
     type(matrix), dimension(:), pointer :: amat,bmat,cmat ! Matrices
     type(matrix_halo), pointer  :: ahalo, chalo ! Haloes
     type(matrix_trans), pointer :: ltrans  ! Points to A/C for type 1/2
     integer, dimension(:), pointer :: bindex ! Index for B for comms
     type(comms_data) :: comms ! Communications type
  end type matrix_mult
!!***

!!****s* matrix_module/trans_remote *
!!  NAME
!!   trans_remote
!!  PURPOSE
!!   Data about remote procs and numbers of pairs to exchange
!!  AUTHOR
!!   D.R.Bowler
!!  SOURCE
!!
  type trans_remote
     ! Maxima
     integer(integ) :: mx_pair
     integer(integ) :: n_rem_node               ! No of halo-nodes
     integer(integ),pointer :: list_rem_node(:) ! Global id of halo-node
     integer(integ),pointer :: nhp_for_node(:)  ! No of halo parts shared 
     integer(integ),pointer :: n_pair(:)        ! No of neig pairs on halo-node
     integer(integ),pointer :: i_pair_addr(:)   ! accumulator for above
     integer(integ),pointer :: i_nd_pair_addr(:)   ! accumulator for above
     integer(integ),pointer :: rem_pair_addr(:) ! For shmems
     integer(integ),pointer :: itran_addr(:)    ! 
  end type trans_remote
!!***

!!****s* matrix_module/pair_data *
!!  NAME
!!   pair_data
!!  PURPOSE
!!   Pointers for transposes - point to an integer array (for
!!   efficient comms)
!!  AUTHOR
!!   D.R.Bowler
!!  SOURCE
!!
  type pair_data
     integer(integ), pointer :: ipart_a(:)  ! mx_iprim*mx_nab
     integer(integ), pointer :: iseq_a(:)   ! mx_iprim*mx_nab
     integer(integ), pointer :: ipart_b(:)  ! mx_iprim*mx_nab
     integer(integ), pointer :: iseq_b(:)   ! mx_iprim*mx_nab
     integer(integ), pointer :: submat(:)
     integer(integ) :: length
  end type pair_data
!!***

contains

!!****f* matrix_module/set_matrix_pointers *
!!
!!  NAME 
!!   set_matrix_pointers
!!  USAGE
!! 
!!  PURPOSE
!!   Establishes pointers for n_nab, i_acc and 
!!   i_seq BEFORE get_naba has been run.  The first two are of length 
!!   number of atoms in partition (and will stay fixed), while the last
!!   can be pointed correctly to the start, and given an upper bound.
!!   The routine set_matrix_pointers2 will establish the correct upper 
!!   bound for i_seq and set the pointer for npxyz
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   24/11/99
!!  MODIFICATION HISTORY
!!   21/06/2001 dave
!!    Added ROBODoc header
!!   08:01, 2003/06/11 dave
!!    Added TM's debugging changes
!!  SOURCE
!!
  subroutine set_matrix_pointers(atoms,mat,ind,part_on_node,mx_part,part_offset)

    ! Passed variables
    integer :: part_on_node, mx_part
    ! Original type(matrix), dimension(part_on_node) :: mat
    type(matrix), dimension(:) :: mat
    integer, dimension(:),target :: ind
    integer, dimension(:) :: atoms, part_offset

    ! Local variables
    integer :: offset, i
    !For debugging
    integer :: isize, isize2

    offset = 0
    do i=1,part_on_node
       offset = part_offset(i)
       if(atoms(i)>mx_part) call cq_abort('set_matrix_pointers: ',i,mx_part)
       if(atoms(i) > 0) then
          isize= atoms(i)
          if(i<part_on_node) then
             isize2 = (part_offset(i+1) - part_offset(i) - 3*isize)/5
          else
             isize2 = (size(ind) - part_offset(i) - 3*isize)/5
          end if
       elseif(atoms(i) == 0) then
          isize=1
          isize2 = 5
       else
          call cq_abort('set_matrix_pointers : size ', atoms(i))
       endif
       mat(i)%n_nab => ind(offset+1:offset+isize)
       offset = offset+isize
       mat(i)%i_acc => ind(offset+1:offset+isize)
       offset = offset+isize
       mat(i)%i_nd_acc => ind(offset+1:offset+isize)
       offset = offset+isize
       mat(i)%i_seq => ind(offset+1:offset+isize2)!mx_part*mat(i)%mx_abs)
       ! Now point at the start of the next data
       ! CHANGED drb 2008/10/06
       !offset = i*(3*mx_part+5*mx_part*mat(i)%mx_abs)
    enddo
    return
  end subroutine set_matrix_pointers
!!***

!!****f* matrix_module/set_matrix_pointers2 *
!!
!!  NAME 
!!   set_matrix_pointers2
!!  USAGE
!! 
!!  PURPOSE
!!   Sets the pointers correctly after get_naba has been called.
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   24/11/99
!!  MODIFICATION HISTORY
!!   21/06/2001 dave
!!    Added ROBODoc header and used cq_abort
!!   08:01, 2003/06/11 dave
!!    Added TM's debugging changes
!!  SOURCE
!!
  subroutine set_matrix_pointers2(atoms,mat,ind,part_on_node,mx_part)

    use GenComms, ONLY: my_barrier

    implicit none

    ! Passed variables
    integer :: part_on_node, mx_part
    ! Original type(matrix), dimension(part_on_node) :: mat
    type(matrix), dimension(:) :: mat
    integer, dimension(:),target :: ind
    integer, dimension(:) :: atoms

    ! Local variables
    integer :: offset, i
    ! For debugging
    integer :: size1, size2, maxlen

!    maxlen = size(ind)
!    if(maxlen<part_on_node*(3*mx_part + 5*mx_part*mat(1)%mx_abs)) &
!         call cq_abort('Overrun: ',maxlen,part_on_node*(3*mx_part + 5*mx_part*mat(1)%mx_abs))
    offset = 0
    do i=1,part_on_node
       !write(io_lun,*) 'Part: ',i
       !call my_barrier
       if(mat(i)%n_atoms>mx_part) then
          call cq_abort('set_matrix_pointers2: too many atoms ',&
               mat(i)%n_atoms,mx_part)
       endif
       if(mat(i)%n_atoms/=atoms(i)) then
          call cq_abort('set_matrix_pointers2: too many atoms(2) ',&
               atoms(i))
       endif
       if(mat(i)%part_nabs>mx_part*mat(i)%mx_abs) then
          call cq_abort('set_matrix_pointers2: too many nabs ',&
               mat(i)%part_nabs,mx_part*mat(i)%mx_abs)
       endif
       mat(i)%array_posn = offset+1
       if(mat(i)%n_atoms >0) then
          size1= mat(i)%n_atoms
       elseif(mat(i)%n_atoms == 0) then
          size1= 1
       else
          call cq_abort('set_matrix_pointers2: n_atoms ',mat(i)%n_atoms)
       endif
       if(mat(i)%part_nabs > 0) then
          size2= mat(i)%part_nabs
       elseif(mat(i)%part_nabs ==0) then
          size2= 1
       else
          call cq_abort('set_matrix_pointers2: part_nabs ',mat(i)%part_nabs)
       endif

       mat(i)%n_nab => ind(offset+1:offset+size1)
       offset = offset+size1
       mat(i)%i_acc => ind(offset+1:offset+size1)
       offset = offset+size1
       mat(i)%i_nd_acc => ind(offset+1:offset+size1)
       offset = offset+size1
       mat(i)%i_seq => ind(offset+1:offset+size2)
       offset = offset+size2
       mat(i)%npxyz => ind(offset+1:offset+3*size2)
       offset = offset+3*size2
       mat(i)%ndimj => ind(offset+1:offset+size2)
       offset = offset+size2
       ! Now point offset to the start of the next partition's data
       ! Actually it already points to it
       !offset = i*(3*mx_part+5*mx_part*mat(i)%mx_abs)
       ! Add an overflow check ?
    enddo
    return
  end subroutine set_matrix_pointers2
!!***

!!****f* matrix_module/set_trans_pointers *
!!
!!  NAME 
!!   set_trans_pointers
!!  USAGE
!! 
!!  PURPOSE
!!   Sets the pointers for the transpose data
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   24/11/99
!!  MODIFICATION HISTORY
!!   21/06/2001 dave
!!    Added ROBODoc header
!!   2008/05/22 ast
!!    Added timers
!!  SOURCE
!!
  subroutine set_trans_pointers(apairind,atrans_rem,apairs)

    ! Module usage
    use datatypes
    use GenComms, ONLY: cq_abort

    ! Passed variables
    integer(integ), target :: apairind(:)
    type(trans_remote) :: atrans_rem
    type(pair_data) :: apairs(:)

    ! Local variables
    integer(integ) :: posn, nr, lenind

    lenind = size(apairind)
    posn = 0
    do nr=1,atrans_rem%n_rem_node
       if(atrans_rem%n_pair(nr) <= 0) call cq_abort('nr, n_pair ',nr, atrans_rem%n_pair(nr))
       apairs(nr)%ipart_a => apairind(posn+1:posn+atrans_rem%n_pair(nr))
       posn = posn+atrans_rem%n_pair(nr)
       apairs(nr)%iseq_a  => apairind(posn+1:posn+atrans_rem%n_pair(nr))
       posn = posn+atrans_rem%n_pair(nr)
       apairs(nr)%ipart_b => apairind(posn+1:posn+atrans_rem%n_pair(nr))
       posn = posn+atrans_rem%n_pair(nr)
       apairs(nr)%iseq_b  => apairind(posn+1:posn+atrans_rem%n_pair(nr))
       posn = posn+atrans_rem%n_pair(nr)
       !write(io_lun,*) 'About to dealloc submat'
       !if(associated(apairs(nr)%submat)) deallocate(apairs(nr)%submat)
       call start_timer(tmr_std_allocation)
       allocate(apairs(nr)%submat(atrans_rem%n_pair(nr)))
       call stop_timer(tmr_std_allocation)
       if(posn>lenind) call cq_abort('Pairs overflow: ',posn,lenind)
    enddo
    !!   ADD AN OVERFLOW CHECK HERE !
    !    nr = atrans_rem%n_rem_node
    !    posn = posn - atrans_rem%n_pair(nr)
    !    if(posn>4*mx_iprim*amat(1)%mx_nab) then
    !      write(io_lun,*) 'Index overflow for pairind ! ', &
    !         posn,4*mx_iprim*amat(1)%mx_nab
    !      stop
    !    endif
    return
  end subroutine set_trans_pointers
!!***

!!****f* matrix_module/allocate_matrix *
!!
!!  NAME 
!!   allocate_matrix
!!  USAGE
!! 
!!  PURPOSE
!!   Allocate memory to matrix derived type
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   24/11/99
!!  MODIFICATION HISTORY
!!   21/06/2001 dave
!!    Added ROBODoc header and cq_abort
!!   2008/05/22 ast
!!    Added timers
!!  SOURCE
!!
  subroutine allocate_matrix(mat,part_on_node,mx_part)

    implicit none

    integer :: mx_part,i,part_on_node,stat
    type(matrix), dimension(:)  :: mat

    call start_timer(tmr_std_allocation)
    do i=1,part_on_node
       allocate(mat(i)%radius(mx_part*mat(i)%mx_abs),STAT=stat)
       if(stat/=0) then
          call cq_abort('alloc_mat: error allocating memory to radius')
       endif
       allocate(mat(i)%i_part(mx_part*mat(i)%mx_abs),STAT=stat)
       if(stat/=0) then
          call cq_abort('alloc_mat: error allocating memory to i_part')
       endif
       !allocate(mat(i)%ndimj(mx_part*mat(i)%mx_abs),STAT=stat)
       !if(stat/=0) then
       !   call cq_abort('alloc_mat: error allocating memory to ndimj')
       !endif
       allocate(mat(i)%onsite(mx_part),STAT=stat)
       if(stat/=0) then
          call cq_abort('alloc_mat: error allocating memory to onsite')
       endif
       allocate(mat(i)%ndimi(mx_part),STAT=stat)
       if(stat/=0) then
          call cq_abort('alloc_mat: error allocating memory to ndimi')
       endif
    enddo
    call stop_timer(tmr_std_allocation)
    return
  end subroutine allocate_matrix
!!***

!!****f* matrix_module/allocate_trans *
!!
!!  NAME 
!!   allocate_trans
!!  USAGE
!! 
!!  PURPOSE
!!   Allocates memory for the transpose derived type
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   24/11/99
!!  MODIFICATION HISTORY
!!   21/06/2001 dave
!!    Added ROBODoc header and cq_abort
!!   2008/05/22 ast
!!    Added timers
!!  SOURCE
!!
  subroutine allocate_trans(trans,mx_iprim)
    
    implicit none

    integer :: mx_iprim,stat
    type(matrix_trans)::trans

    call start_timer(tmr_std_allocation)
    allocate(trans%i_prim(mx_iprim*trans%mx_nab),STAT=stat)
    if(stat/=0) then
       call cq_abort('alloc_trans: error allocating memory to i_prim')
    endif
    allocate(trans%i_beg(trans%mx_halo),trans%i_nd_beg(trans%mx_halo),STAT=stat)
    if(stat/=0) then
       call cq_abort('alloc_trans: error allocating memory to i_beg')
    endif
    allocate(trans%n_hnab(trans%mx_halo),STAT=stat)
    if(stat/=0) then
       call cq_abort('alloc_trans: error allocating memory to i_hnab')
    endif
    call stop_timer(tmr_std_allocation)
    return
  end subroutine allocate_trans
!!***

!!****f* matrix_module/allocate_halo *
!!
!!  NAME 
!!   allocate_halo
!!  USAGE
!! 
!!  PURPOSE
!!   Allocates memory to halo derived type
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   24/11/99
!!  MODIFICATION HISTORY
!!   21/06/2001 dave
!!    Added ROBODoc header and cq_abort
!!   2008/05/22 ast
!!    Added timers
!!  SOURCE
!!
  subroutine allocate_halo(halo,mx_iprim,mx_icover,mx_pcover)

    ! Temp
    use GenComms, ONLY: myid

    implicit none

    integer :: mx_iprim,mx_icover,mx_pcover,stat
    type(matrix_halo):: halo

    call start_timer(tmr_std_allocation)
    allocate(halo%nh_part(halo%mx_part),STAT=stat)
    if(stat/=0) then
       call cq_abort('alloc_halo: error allocating nh_part')
    endif
    allocate(halo%j_beg(halo%mx_part),STAT=stat)
    if(stat/=0) then
       call cq_abort('alloc_halo: error allocating j_beg')
    endif
    allocate(halo%lab_hcell(halo%mx_part),STAT=stat)
    if(stat/=0) then
       call cq_abort('alloc_halo: error allocating lab_hcell')
    endif
    allocate(halo%lab_hcover(halo%mx_part),STAT=stat)
    if(stat/=0) then
       call cq_abort('alloc_halo: error allocating lab_hcover')
    endif
    allocate(halo%j_seq(halo%mx_halo),STAT=stat)
    if(stat/=0) then
       call cq_abort('alloc_halo: error allocating j_seq')
    endif
    allocate(halo%i_h2d(mx_iprim*halo%mx_halo),STAT=stat)
    if(stat/=0) then
       call cq_abort('alloc_halo: error allocating i_h2d')
    endif
    allocate(halo%i_halo(mx_icover),STAT=stat)
    if(stat/=0) then
       call cq_abort('alloc_halo: error allocating i_halo')
    endif
    allocate(halo%i_hbeg(mx_pcover),STAT=stat)
    if(stat/=0) then
       call cq_abort('alloc_halo: error allocating i_hbeg')
    endif
    allocate(halo%ndimj(halo%mx_halo),STAT=stat)
    if(stat/=0) then
       call cq_abort('alloc_halo: error allocating ndimj')
    endif
    allocate(halo%ndimi(mx_iprim),STAT=stat)
    if(stat/=0) then
       call cq_abort('alloc_halo: error allocating ndimi')
    endif
    halo%ndimi = 0
    halo%ndimj = 0
    call stop_timer(tmr_std_allocation)
    return
  end subroutine allocate_halo
!!***

!!****f* matrix_module/allocate_trans_rem *
!!
!!  NAME 
!!   allocate_trans_rem
!!  USAGE
!! 
!!  PURPOSE
!!   Allocates memory to the trans_remote derived type
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   24/11/99
!!  MODIFICATION HISTORY
!!   21/06/2001 dave
!!    Added ROBODoc header and cq_abort
!!   2008/05/22 ast
!!    Added timers
!!  SOURCE
!!
  subroutine allocate_trans_rem(data,mx_neighbour_nodes,size)

    implicit none

    ! Passed variables
    type(trans_remote) :: data
    integer :: mx_neighbour_nodes,size

    ! Local variables
    integer :: stat

    call start_timer(tmr_std_allocation)
    allocate(data%list_rem_node(mx_neighbour_nodes+1),STAT=stat)
    if(stat/=0) then
       call cq_abort('trans_rem: error allocating memory(1)')
    endif
    allocate(data%nhp_for_node(mx_neighbour_nodes+1),STAT=stat)
    if(stat/=0) then
       call cq_abort('trans_rem: error allocating memory(2)')
    endif
    allocate(data%n_pair(mx_neighbour_nodes+1),STAT=stat)
    if(stat/=0) then
       call cq_abort('trans_rem: error allocating memory(3)')
    endif
    allocate(data%i_pair_addr(mx_neighbour_nodes+1),data%i_nd_pair_addr(mx_neighbour_nodes+1),STAT=stat)
    if(stat/=0) then
       call cq_abort('trans_rem: error allocating memory(4)')
    endif
    allocate(data%rem_pair_addr(mx_neighbour_nodes+1),STAT=stat)
    if(stat/=0) then
       call cq_abort('trans_rem: error allocating memory(5)')
    endif
    allocate(data%itran_addr(size),STAT=stat)
    if(stat/=0) then
       call cq_abort('trans_rem: error allocating memory(6)')
    endif
    call stop_timer(tmr_std_allocation)
    return
  end subroutine allocate_trans_rem
!!***

!!****f* matrix_module/deallocate_matrix *
!!
!!  NAME 
!!   deallocate_matrix
!!  USAGE
!! 
!!  PURPOSE
!!   Deallocates the matrix derived type
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   24/11/99
!!  MODIFICATION HISTORY
!!   21/06/2001 dave
!!    Added ROBODoc header and cq_abort
!!   2008/05/22 ast
!!    Added timers
!!  SOURCE
!!
  subroutine deallocate_matrix(mat,part_on_node)

    implicit none

    integer :: i,part_on_node,stat
    type(matrix), dimension(:) :: mat

    call start_timer(tmr_std_allocation)
    do i=part_on_node,1,-1
       deallocate(mat(i)%ndimi,STAT=stat)
       if(stat/=0) then
          call cq_abort('dealloc_matrix: error deallocating memory(1a)')
       endif
       deallocate(mat(i)%onsite,STAT=stat)
       if(stat/=0) then
          call cq_abort('dealloc_matrix: error deallocating memory(1)')
       endif
       !deallocate(mat(i)%ndimj,STAT=stat)
       !if(stat/=0) then
       !   call cq_abort('dealloc_matrix: error deallocating memory(2a)')
       !endif
       deallocate(mat(i)%i_part,STAT=stat)
       if(stat/=0) then
          call cq_abort('dealloc_matrix: error deallocating memory(2)')
       endif
       deallocate(mat(i)%radius,STAT=stat)
       if(stat/=0) then
          call cq_abort('dealloc_matrix: error deallocating memory(3)')
       endif
    enddo
    call stop_timer(tmr_std_allocation)
    return
  end subroutine deallocate_matrix
!!***

!!****f* matrix_module/deallocate_trans *
!!
!!  NAME 
!!   deallocate_trans
!!  USAGE
!! 
!!  PURPOSE
!!   Deallocates the trans derived type
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   24/11/99
!!  MODIFICATION HISTORY
!!   21/06/2001 dave
!!    Added ROBODoc header and cq_abort
!!   2008/05/22 ast
!!    Added timers
!!  SOURCE
!!
  subroutine deallocate_trans(trans)

    implicit none

    type(matrix_trans):: trans
    integer :: stat

    call start_timer(tmr_std_allocation)
    deallocate(trans%n_hnab,STAT=stat)
    if(stat/=0) then
       call cq_abort('dealloc_trans: error deallocating memory(1)')
    endif
    deallocate(trans%i_nd_beg,trans%i_beg,STAT=stat)
    if(stat/=0) then
       call cq_abort('dealloc_trans: error deallocating memory(2)')
    endif
    deallocate(trans%i_prim,STAT=stat)
    if(stat/=0) then
       call cq_abort('dealloc_trans: error deallocating memory(3)')
    endif
    call stop_timer(tmr_std_allocation)
    return
  end subroutine deallocate_trans
!!***

!!****f* matrix_module/deallocate_halo *
!!
!!  NAME 
!!   deallocate_halo
!!  USAGE
!! 
!!  PURPOSE
!!   deallocates halo derived type
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   24/11/99
!!  MODIFICATION HISTORY
!!   21/06/2001 dave
!!    Added ROBODoc header and cq_abort
!!    Bug fix - added missing bracket
!!   2008/05/22 ast
!!    Added timers
!!  SOURCE
!!
  subroutine deallocate_halo(halo)

    implicit none

    type(matrix_halo):: halo
    integer :: stat

    stat = 0
    !if(associated(halo%ndimi)) deallocate(halo%ndimi,STAT=stat)
    !if(stat/=0) call cq_abort('dealloc_halo: error deallocating memory from ndimi',stat)
    !if(associated(halo%ndimj)) deallocate(halo%ndimj,STAT=stat)
    !if(stat/=0) call cq_abort('dealloc_halo: error deallocating memory from ndimj',stat)
    !if(associated(halo%i_hbeg)) deallocate(halo%i_hbeg,STAT=stat)
    !if(stat/=0) call cq_abort('dealloc_halo: error deallocating memory from i_hbeg',stat)
    !if(associated(halo%i_halo)) deallocate(halo%i_halo,STAT=stat)
    !if(stat/=0) call cq_abort('dealloc_halo: error deallocating memory(1)')
    !if(associated(halo%i_h2d)) deallocate(halo%i_h2d,STAT=stat)
    !if(stat/=0) call cq_abort('dealloc_halo: error deallocating memory(2)')
    !if(associated(halo%j_seq)) deallocate(halo%j_seq,STAT=stat)
    !if(stat/=0) call cq_abort('dealloc_halo: error deallocating memory(3)')
    !if(associated(halo%lab_hcover)) deallocate(halo%lab_hcover,STAT=stat)
    !if(stat/=0) call cq_abort('dealloc_halo: error deallocating memory(4)')
    !if(associated(halo%lab_hcell)) deallocate(halo%lab_hcell,STAT=stat)
    !if(stat/=0) call cq_abort('dealloc_halo: error deallocating memory(5)')
    !if(associated(halo%j_beg)) deallocate(halo%j_beg,STAT=stat)
    !if(stat/=0) call cq_abort('dealloc_halo: error deallocating memory(6)')
    !if(associated(halo%nh_part)) deallocate(halo%nh_part,STAT=stat)
    !if(stat/=0) call cq_abort('dealloc_halo: error deallocating memory(7)')
    call start_timer(tmr_std_allocation)
    deallocate(halo%ndimi,STAT=stat)
    if(stat/=0) call cq_abort('dealloc_halo: error deallocating memory from ndimi',stat)
    deallocate(halo%ndimj,STAT=stat)
    if(stat/=0) call cq_abort('dealloc_halo: error deallocating memory from ndimj',stat)
    deallocate(halo%i_hbeg,STAT=stat)
    if(stat/=0) call cq_abort('dealloc_halo: error deallocating memory from i_hbeg',stat)
    deallocate(halo%i_halo,STAT=stat)
    if(stat/=0) call cq_abort('dealloc_halo: error deallocating memory(1)')
    deallocate(halo%i_h2d,STAT=stat)
    if(stat/=0) call cq_abort('dealloc_halo: error deallocating memory(2)')
    deallocate(halo%j_seq,STAT=stat)
    if(stat/=0) call cq_abort('dealloc_halo: error deallocating memory(3)')
    deallocate(halo%lab_hcover,STAT=stat)
    if(stat/=0) call cq_abort('dealloc_halo: error deallocating memory(4)')
    deallocate(halo%lab_hcell,STAT=stat)
    if(stat/=0) call cq_abort('dealloc_halo: error deallocating memory(5)')
    deallocate(halo%j_beg,STAT=stat)
    if(stat/=0) call cq_abort('dealloc_halo: error deallocating memory(6)')
    deallocate(halo%nh_part,STAT=stat)
    if(stat/=0) call cq_abort('dealloc_halo: error deallocating memory(7)')
    call stop_timer(tmr_std_allocation)
    return
  end subroutine deallocate_halo
!!***

!!****f* matrix_module/deallocate_trans_rem *
!!
!!  NAME 
!!   deallocate_trans_rem
!!  USAGE
!! 
!!  PURPOSE
!!   deallocates trans_remote derived type
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   24/11/99
!!  MODIFICATION HISTORY
!!   21/06/2001 dave
!!    Added ROBODoc header and cq_abort
!!   2008/05/22 ast
!!    Added timers
!!  SOURCE
!!
  subroutine deallocate_trans_rem(data)

    ! Passed variables
    type(trans_remote) :: data

    ! Local variables
    integer :: stat

    call start_timer(tmr_std_allocation)
    deallocate(data%itran_addr,data%rem_pair_addr,data%i_pair_addr,data%i_nd_pair_addr, &
         data%n_pair,data%nhp_for_node,data%list_rem_node,STAT=stat)
    if(stat/=0) then
       call cq_abort('dealloc_trans_rem: error deallocating')
    endif
    call stop_timer(tmr_std_allocation)

    return
  end subroutine deallocate_trans_rem
!!***

!!****f* matrix_module/allocate_comms_data *
!!
!!  NAME 
!!   allocate_comms_data
!!  USAGE
!! 
!!  PURPOSE
!!   Allocates the comms_data derived type
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   24/11/99
!!  MODIFICATION HISTORY
!!   21/06/2001 dave
!!    Added ROBODoc header and cq_abort
!!   2008/05/22 ast
!!    Added timers
!!  SOURCE
!!
  subroutine allocate_comms_data(a_b_c,mx_part,mx_nab_nodes,mx_nponn)

    implicit none

    ! Passed variables
    integer :: mx_part,mx_nab_nodes,mx_nponn
    type(comms_data) :: a_b_c

    ! Local variables
    integer :: stat

    call start_timer(tmr_std_allocation)
    allocate(a_b_c%np_send(mx_nab_nodes),STAT=stat) 
    if(stat/=0) then
       call cq_abort('alloc_comms_data: error allocating memory(1)')
    endif
    allocate(a_b_c%pl_send(mx_part,mx_nab_nodes),STAT=stat)
    if(stat/=0) then
       call cq_abort('alloc_comms_data: error allocating memory(2)')
    endif
    allocate(a_b_c%neigh_node_list(mx_part),STAT=stat)
    if(stat/=0) then
       call cq_abort('alloc_comms_data: error allocating memory(3)')
    endif
    allocate(a_b_c%ilen2rec(mx_nponn,mx_nab_nodes),a_b_c%ilen3rec(mx_nponn,mx_nab_nodes),STAT=stat) 
    if(stat/=0) then
       call cq_abort('alloc_comms_data: error allocating memory(4)')
    endif
    allocate(a_b_c%istart(mx_nponn,mx_nab_nodes),STAT=stat)
    if(stat/=0) then
       call cq_abort('alloc_comms_data: error allocating memory(5)')
    endif
    allocate(a_b_c%ncomm(mx_nab_nodes),STAT=stat) 
    if(stat/=0) then
       call cq_abort('alloc_comms_data: error allocating memory(6)')
    endif
    call stop_timer(tmr_std_allocation)
  end subroutine allocate_comms_data
!!***

!!****f* matrix_module/deallocate_comms_data *
!!
!!  NAME 
!!   deallocate_comms_data
!!  USAGE
!! 
!!  PURPOSE
!!   deallocates the comms_data derived type
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   24/11/99
!!  MODIFICATION HISTORY
!!   21/06/2001 dave
!!    Added ROBODoc header and cq_abort
!!  SOURCE
!!
  subroutine deallocate_comms_data(a_b_c)

    implicit none

    type(comms_data) :: a_b_c

    integer :: stat

    
    stat = 0
    deallocate(a_b_c%ncomm,STAT=stat)
    if(stat/=0) call cq_abort('dealloc_comms_data: error 1 deallocating memory')
    deallocate(a_b_c%istart,STAT=stat)
    if(stat/=0) call cq_abort('dealloc_comms_data: error 2 deallocating memory')
    deallocate(a_b_c%ilen2rec,a_b_c%ilen3rec,STAT=stat)
    if(stat/=0) call cq_abort('dealloc_comms_data: error 3 deallocating memory')
    deallocate(a_b_c%neigh_node_list,STAT=stat)
    if(stat/=0) call cq_abort('dealloc_comms_data: error 4 deallocating memory')
    deallocate(a_b_c%pl_send,STAT=stat)
    if(stat/=0) call cq_abort('dealloc_comms_data: error 5 deallocating memory')
    deallocate(a_b_c%np_send,STAT=stat)
    if(stat/=0) call cq_abort('dealloc_comms_data: error 6 deallocating memory')
    !if(associated(a_b_c%ncomm)) deallocate(a_b_c%ncomm,STAT=stat)
    !if(stat/=0) call cq_abort('dealloc_comms_data: error 1 deallocating memory')
    !if(associated(a_b_c%istart)) deallocate(a_b_c%istart,STAT=stat)
    !if(stat/=0) call cq_abort('dealloc_comms_data: error 2 deallocating memory')
    !if(associated(a_b_c%ilen2rec)) deallocate(a_b_c%ilen2rec,a_b_c%ilen3rec,STAT=stat)
    !if(stat/=0) call cq_abort('dealloc_comms_data: error 3 deallocating memory')
    !if(associated(a_b_c%neigh_node_list)) deallocate(a_b_c%neigh_node_list,STAT=stat)
    !if(stat/=0) call cq_abort('dealloc_comms_data: error 4 deallocating memory')
    !if(associated(a_b_c%pl_send)) deallocate(a_b_c%pl_send,STAT=stat)
    !if(stat/=0) call cq_abort('dealloc_comms_data: error 5 deallocating memory')
    !if(associated(a_b_c%np_send)) deallocate(a_b_c%np_send,STAT=stat)
    !if(stat/=0) call cq_abort('dealloc_comms_data: error 6 deallocating memory')
  end subroutine deallocate_comms_data
!!***
end module matrix_module
