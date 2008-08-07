! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module atoms
! ------------------------------------------------------------------------------
! Code area 8: Indexing
! ------------------------------------------------------------------------------

!!****h* Conquest/atoms *
!!  NAME
!!   atoms
!!  PURPOSE
!!   Stores atom positions, distributions among processors
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   14/05/01
!!  MODIFICATION HISTORY
!!   18/03/2002 dave
!!    Added RCS Id and Log tags and static tag for object file id
!!   2004/11/10 drb
!!    Changed dimensions to use global and mx_, and added allocation call
!!   2008/02/04 08:24 dave
!!    Changed for output to file not stdout
!!   2008/05/16 ast
!!    Added some timers
!!  SOURCE
!!
module atoms

  use datatypes
  use global_module, ONLY: io_lun
  use timer_stdclocks_module, ONLY: start_timer,stop_timer,start_timer, stop_timer, tmr_std_indexing,tmr_std_allocation

  implicit none
  save

  ! RCS tag for object file identification 
  character(len=80), private :: RCSid = "$Id$"

  integer :: n_my_atoms

  integer, dimension(:), allocatable :: index_my_atoms, node_doing_atom, atom_number_on_node, n_atoms_on_node
  integer, dimension(:,:), allocatable :: atoms_on_node


!!***

contains

! -----------------------------------------------------------
! Subroutine distribute_atoms
! -----------------------------------------------------------

!!****f* atoms/distribute_atoms *
!!
!!  NAME 
!!   distribute_atoms
!!  USAGE
!! 
!!  PURPOSE
!!   Identifies which node is dealing with which atom and 
!!   creates tables for switch back and forth between the
!!   node-list and global list.
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   C.M.Goringe/D.R.Bowler
!!  CREATION DATE
!!   23/02/95
!!  MODIFICATION HISTORY
!!   12/9/95 by Chris Goringe. Atoms no longer 
!!     distributed by domain. 
!!   22/3/96 by Mike Gillan. Atoms now assigned 
!!     to nodes in sequential groups: for any 
!!     two atoms i1 and i2 on nodes n1 and n2 resp., 
!!     then i2 .gt. i1 implies n2 .ge. n1. 
!!   02/05/00 by Dave Bowler.  Since atoms are now distributed 
!!     by partitions, this routine must reflect that !
!!   25/05/2001 dave
!!    Added ROBODoc header, stripped subroutine call, placed
!!    inside atoms.module
!!   11:23, 12/11/2004 dave 
!!    Added allocation call
!!   2008/05/16 ast
!!    Added timer
!!  SOURCE
!!
  subroutine distribute_atoms(inode, ionode) 
    
    use global_module, ONLY: id_glob, numprocs, iprint_init
    use group_module, ONLY: parts
    
    implicit none 

    ! Shared variables 
    integer :: inode, ionode

    ! Local variables 
    integer :: n_per_node, n_first, ia, node_doing, i, ni, np
    integer :: ind_part


    call start_timer(tmr_std_indexing)

    if(iprint_init>1.AND.inode==ionode) write(io_lun,fmt='(10x,"Allocating memory for distribute_atom")')
    call allocate_distribute_atom
    ia = 0
    do node_doing = 1,numprocs
       ! Start by working out how many atoms each node owns
       n_atoms_on_node(node_doing) = 0
       do np = 1,parts%ng_on_node(node_doing)
          ind_part = parts%ngnode(parts%inode_beg(node_doing)+np-1)
          n_atoms_on_node(node_doing) = n_atoms_on_node(node_doing) + &
               parts%nm_group(ind_part)
       enddo
       i=0
       ! Now distribute them
       do np = 1,parts%ng_on_node(node_doing)
          ind_part = parts%ngnode(parts%inode_beg(node_doing)+np-1)
          if(parts%nm_group(ind_part) > 0) then
             do ni = 1,parts%nm_group(ind_part)
                ia = id_glob(parts%icell_beg(ind_part)+ni-1)
                i = i+1      
                atoms_on_node(i, node_doing) = ia ! global id
                atom_number_on_node(ia) = i       ! primary number
                node_doing_atom(ia) = node_doing  ! locator
             enddo
          end if
       enddo
    enddo
    n_my_atoms = n_atoms_on_node(INODE)
    do i = 1, n_my_atoms 
       index_my_atoms(i) = atoms_on_node(i, INODE) 
    enddo
    call stop_timer(tmr_std_indexing)
    return 
  end subroutine distribute_atoms
!!*** 

! ------------------------------------------------------------------------------
! Subroutine allocate_distribute_atom
! ------------------------------------------------------------------------------

!!****f* atoms/allocate_distribute_atom *
!!
!!  NAME 
!!   allocate_distribute_atom
!!  USAGE
!! 
!!  PURPOSE
!!   Allocates dynamic variables associated with atom distribution across 
!!   processors (index_my_atoms and atoms_on_node)
!!  INPUTS
!! 
!! 
!!  USES
!!   global_module, primary_module
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!!   mid-2000 ?
!!  MODIFICATION HISTORY
!!   22/03/2002 dave
!!    Added ROBODoc header
!!   11:23, 12/11/2004 dave 
!!    Removed from distribute_atom
!!   2008/05/15 ast
!!    Added timer
!!  SOURCE
!!
  subroutine allocate_distribute_atom

    use global_module, ONLY: numprocs, ni_in_cell, area_index
    use primary_module, ONLY: bundle
    use GenComms, ONLY: cq_abort
    use memory_module, ONLY: reg_alloc_mem, type_int
    
    implicit none

    ! Local variables
    integer :: stat

    call start_timer(tmr_std_allocation)

    allocate(index_my_atoms(bundle%mx_iprim),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to index_my_atoms !')
    call reg_alloc_mem(area_index, bundle%mx_iprim, type_int)
    allocate(atoms_on_node(bundle%mx_iprim,numprocs),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to atoms_on_node !')
    call reg_alloc_mem(area_index, bundle%mx_iprim, type_int)
    allocate(node_doing_atom(ni_in_cell), atom_number_on_node(ni_in_cell), n_atoms_on_node(numprocs), STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to node_doing_atom !')
    call reg_alloc_mem(area_index, 2*ni_in_cell+numprocs, type_int)
    call stop_timer(tmr_std_allocation)
    return
  end subroutine allocate_distribute_atom
!!***
end module atoms
