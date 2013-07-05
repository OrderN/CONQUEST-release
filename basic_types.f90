! --------------------------------------------------------------------
! Module basic_types
! --------------------------------------------------------------------

!!****h* Conquest/basic_types *
!!  NAME
!!   basic_types - defines derived types
!!  PURPOSE
!!   Define derived types for small groups, primary sets and covering
!!     sets for atoms and grid points
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   20/04/00
!!  MODIFICATION HISTORY
!! 
!!***
module basic_types

  ! Module usage
  use datatypes

  implicit none

  ! Group refers to generic small organisational group (e.g. block or 
  ! partition) 
  ! Member refers to either atoms or integration grid points

!!****s* basic_types/group_set *
!!  NAME
!!   group_set -- small group type
!!  PURPOSE
!!   Defines a derived type to contain all variables required for small
!!     groups, including different indexing types
!!  AUTHOR
!!   D.R.Bowler
!!  SOURCE
!!
  type group_set
    ! Maxima 
    integer :: mx_ngonn,mx_gcell  ! Maximum groups on a node and in cell
    integer :: mx_gedge   ! Max groups on sim cell edge
    integer :: mx_mem_grp ! Max members in group (was mx_part)
    ! Scalars
    integer :: ngcellx,ngcelly,ngcellz  ! Number of groups in cell sides
    ! Arrays
    integer, pointer :: ng_on_node(:)   ! Number of groups on node
    integer, pointer :: inode_beg(:)    ! Where the groups for a node start
    integer, pointer :: ngnode(:)       ! CC label of group (node order)
    integer, pointer :: i_cc2node(:)    ! Node owning group (CC label)
    integer, pointer :: i_cc2seq(:)     ! Seq. no of group on node (CC label)
    integer, pointer :: nm_group(:)     ! Number of members in a group (CC)
    integer, pointer :: icell_beg(:)    ! Start of member data for a group (CC)
    integer, pointer :: inv_ngnode(:)   ! Inverse of ngnode
  end type group_set
!!***

!!****s* basic_types/primary_set *
!!  NAME
!!   primary_set -- primary set type
!!  PURPOSE
!!   Defines a derived type that includes all information about 
!!     the primary set made up of small groups
!!  AUTHOR
!!   D.R.Bowler
!!  SOURCE
!!
  type primary_set  ! Contains information about a primary set
    ! Define the kind of groups I'm made from 
    type(group_set), pointer :: groups
    ! Scalars
    integer :: mx_iprim ! Maximum members of a primary set
    integer :: n_prim  ! Number of members in primary set
    integer :: mx_ngonn  ! Maximum groups in a primary set
    integer :: groups_on_node ! Actual groups in a primary set
    integer :: nx_origin,ny_origin,nz_origin ! Origin group of primary set
    integer :: nw_primx,nw_primy,nw_primz  ! Widths of primary set
    integer :: nleftx,nlefty,nleftz        ! Left-spans of primary set
    ! Properties of the groups (dim mx_ngonn)
    integer, pointer :: idisp_primx(:) ! Displacement of group from origin  
    integer, pointer :: idisp_primy(:)
    integer, pointer :: idisp_primz(:)
    integer, pointer :: nm_nodgroup(:)  ! Members in a group
    integer, pointer :: nm_nodbeg(:)    ! Where the members for a group start
    ! Properties of the members (dim mx_iprim)
    integer, pointer :: iprim_seq(:)
    integer, pointer :: ig_prim(:)
    real(double), pointer :: xprim(:)  ! Location of members
    real(double), pointer :: yprim(:)
    real(double), pointer :: zprim(:)
    integer, pointer :: species(:)     ! Species of members
  end type primary_set
!!***

!!****s* basic_types/cover_set
!!  NAME
!!   cover_set -- covering set type
!!  PURPOSE
!!   Defines all variables associated with a covering set of one small
!!     group type around the primary set of (possibly another) small
!!     group type
!!  AUTHOR
!!   D.R.Bowler
!!  MODIFICATION HISTORY
!!   2013/07/01 M.Arita
!!    Added ig_cover along with MD implementation
!!  SOURCE
!!
  type cover_set
    ! Define the kind of groups I'm made from 
    type(group_set), pointer :: groups
    ! And the primary set around which I'm built
    type(primary_set), pointer :: prim
    ! Cutoff
    real(double) :: rcut
    ! Maxima
    integer :: mx_gcover                ! max no. of groups in covering set
    integer :: mx_mcover                ! max no. of members in covering set
    ! Scalars
    integer :: ng_cover                 ! no. of groups in CS
    integer :: ncoverx,ncovery,ncoverz  ! x-, y- and z-widths of CS
    integer :: nspanlx,nspanly,nspanlz  ! x-, y- and z-left spans of CS
    integer :: nx_origin,ny_origin,nz_origin ! Origin group of cover set
    ! Arrays (first mx_pcover, then mx_icover)
    integer, pointer :: lab_cell(:)     ! sim-cell gp (CC) equiv to CS gp (PG)
    integer, pointer :: lab_cover(:)    ! CC label of CS group (PG)
    ! TM temporary ????
    integer, pointer :: inv_lab_cover(:)    ! inverse func. of lab_cover CC -> PG
    integer, pointer :: n_ing_cover(:)  ! no. of members in CS partition i (PG)
    integer, pointer :: icover_ibeg(:)  ! accumulator for n_ing_cover
    real(double), pointer :: xcover(:)  ! coords of member (PG order) in CS
    real(double), pointer :: ycover(:)  ! 
    real(double), pointer :: zcover(:)  ! 
    integer, pointer :: ncover_rem(:)   ! Values of ncover from ALL remote nodes
    integer, pointer :: spec_cover(:)   ! List of species of atoms in CS
    integer, pointer :: iprim_group(:)  ! Gives GCS label for a primary member
    integer, pointer :: ig_cover(:)     ! CS --> global label
  end type cover_set
!!***
end module basic_types

