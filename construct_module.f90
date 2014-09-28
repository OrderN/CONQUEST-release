!!****h* Conquest/construct_module
!!  NAME
!!   construct_module
!!  PURPOSE
!!   Aims to provide subroutines that act a little like
!!   constructors for derived types - i.e. allocates memory
!!   for them and initialises one or two variables. Still in
!!   a rudimentary form.
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   26/04/00
!!  MODIFICATION HISTORY
!!   28/05/2001 dave
!!    ROBODoc header, indenting
!!***
module construct_module
  
  !use output_module, only: sub_enter_output, sub_leave_output
  use timer_module,  only: start_timer, stop_timer, cq_timer

contains
!!****f* construct_module/init_group *
!!
!!  NAME 
!!   init_group - initialises groups of atoms/grid points
!!  USAGE
!! 
!!  PURPOSE
!!   Initialise the derived type associated with groups of objects - 
!!   typically either blocks of grid points or partitions of atoms.
!!   This group_set type (defined in basic_types) has global information
!!   stored in it - numbers of groups in cell side, processor ownership of
!!   groups etc.
!!  INPUTS
!!   type(group_set) :: groups - the derived type
!!   integer :: mx_ngonn,mx_gedge,mx_gcell,mx_mem_grp,mx_node
!!    Maxima - groups on a node, groups in an edge, groups in the cell, 
!!    members in a group, nodes
!!  USES
!!   basic_types, group_module
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   26/04/00
!!  MODIFICATION HISTORY
!!   28/05/2001 dave
!!    ROBODoc header and better documentation
!!  SOURCE
!!
  subroutine init_group(groups,mx_ngonn,mx_gedge,mx_gcell,mx_mem_grp,mx_node)
    
    ! Module usage
    use basic_types
    use group_module, only: allocate_group_set
    
    implicit none
    
    ! Passed variables
    type(group_set) :: groups
    integer         :: mx_ngonn,mx_gedge,mx_gcell,mx_mem_grp,mx_node

    ! Local variables
    type(cq_timer)  :: tmr_std_loc

!****lat<$
    call start_timer(t=tmr_std_loc,who='init_group',where=1,level=4)
!****lat>$
 
    groups%mx_ngonn=mx_ngonn
    groups%mx_gedge=mx_gedge
    groups%mx_gcell=mx_gcell
    groups%mx_mem_grp=mx_mem_grp
    call allocate_group_set(groups,mx_node)

!****lat<$
    call stop_timer(t=tmr_std_loc,who='init_group')
!****lat>$

    return
  end subroutine init_group
!!***
  
!!****f* construct_module/init_primary *
!!
!!  NAME 
!!   init_primary - initialises a primary set
!!  USAGE
!! 
!!  PURPOSE
!!   Initialise the derived type associated with a primary set - that is
!!   the set of groups associated with a given processor.  This primary_set
!!   type (defined in basic_types) is purely local, and refers only to the
!!   processor's groups (a bundle of partitions or a domain of blocks).
!!  INPUTS
!!   type(primary_set) :: prim - the derived type
!!   integer :: mx_iprim, mx_ngonn - maxima (members in set, groups on a node)
!!   logical :: members - do we store information about the members of the groups
!!  USES
!!   basic_types, group_module
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   26/04/00
!!  MODIFICATION HISTORY
!!   28/05/2001 dave
!!    ROBODoc header
!!  SOURCE
!!
  subroutine init_primary(prim,mx_iprim,mx_ngonn,members)
    
    ! Module usage
    use basic_types
    use primary_module, only: allocate_primary_set
    
    implicit none
    
    ! Passed variables
    type(primary_set) :: prim
    integer           :: mx_iprim, mx_ngonn
    logical           :: members

    ! Local variables
    type(cq_timer)    :: tmr_std_loc

!****lat<$
    call start_timer(t=tmr_std_loc,who='init_primary',where=1,level=2)
!****lat>$

    prim%mx_iprim = mx_iprim
    prim%mx_ngonn = mx_ngonn
    call allocate_primary_set(prim,members)

!****lat<$
    call stop_timer(t=tmr_std_loc,who='init_primary')
!****lat>$

    return
  end subroutine init_primary
!!***
  
!!****f* construct_module/init_cover *
!!
!!  NAME 
!!   init_cover
!!  USAGE
!! 
!!  PURPOSE
!!   Initialise the derived type associated with a covering set - that is
!!   an orthorhombic set of groups which encompasses all groups within a
!!   specified radius of a primary set.  This cover_set type (defined in 
!!   basic_types) is essentially local, provided the unit cell is big 
!!   enough.
!!  INPUTS
!!   type(cover_set) :: set - the derived type
!!   integer :: mx_mcover, mx_gcover, mx_node (maxima - members in covering
!!    set, groups in cover set, nodes)
!!   logical :: members - do we find information about the members in set
!!  USES
!!   basic_types, cover_module
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   26/04/00
!!  MODIFICATION HISTORY
!!   28/05/2001 dave
!!    ROBODoc header, more documentation
!!  SOURCE
!!
  subroutine init_cover(set,mx_mcover,mx_gcover,mx_node,members)
    
    ! Module usage
    use basic_types
    use cover_module, only: allocate_cs
    
    implicit none
    
    ! Passed variables
    type(cover_set) :: set
    integer         :: mx_mcover,mx_gcover,mx_node
    logical         :: members
    
    ! Local variables
    type(cq_timer)   :: tmr_std_loc


!****lat<$
    call start_timer(t=tmr_std_loc,who='init_cover',where=1,level=4)
!****lat>$

    set%mx_mcover = mx_mcover
    set%mx_gcover = mx_gcover
    call allocate_cs(set,mx_node,members)

!****lat<$
    call stop_timer(t=tmr_std_loc,who='init_cover')
!****lat>$

    return
  end subroutine init_cover
!!***
    
end module construct_module
