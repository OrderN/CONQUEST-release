! $Id$
! -----------------------------------------------------------
! Module group_module
! -----------------------------------------------------------
! Code area 8: indexing
! -----------------------------------------------------------

!!****h* Conquest/group_module
!!  NAME
!!   group_module
!!  PURPOSE
!!   This deals with, and contains variables concerning, the groups 
!!   in the unit cell (i.e. partitions and blocks of integration 
!!   grid points)
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   20/04/00
!!  MODIFICATION HISTORY
!!   30/05/2001 dave
!!    Added ROBODoc headers
!!   20/06/2001 dave
!!    Added RCS Id and Log tags and used cq_abort throughout
!!   2008/05/16 ast
!!    Added timers
!!  SOURCE
!!
module group_module

  ! Module usage
  use datatypes
  use basic_types
  use timer_stdclocks_module, ONLY: start_timer,stop_timer,tmr_std_indexing,tmr_std_allocation

  implicit none
  save

  type(group_set) :: parts  ! Partitions of atoms
  type(group_set) :: blocks ! Blocks of integration grid points

  integer :: part_method
  integer, parameter :: HILBERT = 1
  integer, parameter :: PYTHON = 2
!!***

contains

!!****f* group_module/make_cc2 *
!!
!!  NAME 
!!   make_cc2
!!  USAGE
!! 
!!  PURPOSE
!!   This subroutine constructs arrays which convert from absolute number
!!   of a group to it's "sequence" number on it's node (and tells you 
!!   which node is responsible for it)
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   20/04/00
!!  MODIFICATION HISTORY
!!   30/05/2001 dave
!!    Added ROBODoc header
!!   2008/05/16 ast
!!    Added timer
!!  SOURCE
!!
  subroutine make_cc2(groups,numprocs)

    ! Module usage
    use datatypes
    use basic_types

    implicit none

    ! Passed variables
    type(group_set) :: groups
    integer :: numprocs

    ! Local variables
    integer :: nnd,np,ind_group

    call start_timer(tmr_std_indexing)
    do nnd=1,numprocs  ! Loop over processors
      if(groups%ng_on_node(nnd).gt.0) then  
        do np=1,groups%ng_on_node(nnd)  ! Loop over groups on the node
          ind_group=groups%ngnode(groups%inode_beg(nnd)+np-1)  
          groups%i_cc2node(ind_group)=nnd
          groups%i_cc2seq(ind_group)=np
        enddo
      endif
    enddo   
    call stop_timer(tmr_std_indexing)
  end subroutine make_cc2
!!***

!!****f* group_module/allocate_group_set *
!!
!!  NAME 
!!   allocate_group_set
!!  USAGE
!! 
!!  PURPOSE
!!   Allocates memory and assigns maxima to the group_set type
!!  INPUTS
!! 
!! 
!!  USES
!!   datatypes, basic_types
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   20/04/00
!!  MODIFICATION HISTORY
!!   30/05/2001 dave
!!    Added ROBODoc header
!!   2008/05/16 ast
!!    Added timer
!!   2009/11/03 16:43 dave
!!    Added memory registration
!!  SOURCE
!!
  subroutine allocate_group_set(groups,mx_node)

    ! Module usage
    use datatypes
    use basic_types
    use GenComms, ONLY: cq_abort
    use memory_module, ONLY: reg_alloc_mem, type_int
    use global_module, ONLY: area_index

    implicit none

    ! Passed variables
    type(group_set) :: groups
    integer :: mx_node

    ! Local variables
    integer :: stat

    call start_timer(tmr_std_allocation)
    allocate(groups%ng_on_node(mx_node),STAT=stat)
    if(stat/=0) then
       call cq_abort('alloc_gp: error allocating memory to ng_on_node !')
    endif
    allocate(groups%inode_beg(mx_node),STAT=stat)
    if(stat/=0) then
       call cq_abort('alloc_gp: error allocating memory to inode_beg !')
    endif
    call reg_alloc_mem(area_index,2*mx_node,type_int)
    allocate(groups%ngnode(groups%mx_gcell),STAT=stat)
    if(stat/=0) then
       call cq_abort('alloc_gp: error allocating memory to ngnode !')
    endif
    allocate(groups%i_cc2node(groups%mx_gcell),STAT=stat)
    if(stat/=0) then
       call cq_abort('alloc_gp: error allocating memory to i_cc2node !')
    endif
    allocate(groups%i_cc2seq(groups%mx_gcell),STAT=stat)
    if(stat/=0) then
       call cq_abort('alloc_gp: error allocating memory to i_cc2seq !')
    endif
    allocate(groups%nm_group(groups%mx_gcell),STAT=stat)
    if(stat/=0) then
       call cq_abort('alloc_gp: error allocating memory to nm_group !')
    endif
    allocate(groups%icell_beg(groups%mx_gcell),STAT=stat)
    if(stat/=0) then
       call cq_abort('alloc_gp: error allocating memory to icell_beg !')
    endif
    allocate(groups%inv_ngnode(groups%mx_gcell),STAT=stat)
    if(stat/=0) then
       call cq_abort('alloc_gp: error allocating memory to inv_ngnode !')
    endif
    call reg_alloc_mem(area_index,6*groups%mx_gcell,type_int)
    call stop_timer(tmr_std_allocation)
    return
  end subroutine allocate_group_set
!!***

!!****f* group_module/deallocate_group_set *
!!
!!  NAME 
!!   deallocate_group_set
!!  USAGE
!! 
!!  PURPOSE
!!   Deallocates memory associated with the group_set type
!!  INPUTS
!! 
!! 
!!  USES
!!   datatypes, basic_types
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   20/04/00
!!  MODIFICATION HISTORY
!!   30/05/2001 dave
!!    Added ROBODoc header
!!   2008/05/16 ast
!!    Added timer
!!   2009/11/03 16:43 dave
!!    Added memory registration
!!  SOURCE
!!
  subroutine deallocate_group_set(groups)

    use basic_types
    use GenComms, ONLY: cq_abort
    use memory_module, ONLY: reg_dealloc_mem, type_int
    use global_module, ONLY: area_index

    implicit none

    ! Passed variables
    type(group_set) :: groups

    ! Local variables
    integer :: stat

    call start_timer(tmr_std_allocation)
    deallocate(groups%inv_ngnode,groups%icell_beg,groups%i_cc2seq, &
         groups%i_cc2node,groups%ngnode,groups%inode_beg, &
         groups%ng_on_node,STAT=stat)
    if(stat/=0) then
       call cq_abort('dealloc_gp: error deallocating group_set !')
    endif
    call reg_dealloc_mem(area_index,6*groups%mx_gcell,type_int)
    call stop_timer(tmr_std_allocation)
    return
  end subroutine deallocate_group_set
!!***
end module group_module
