! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module bucket_module
! ------------------------------------------------------------------------------
! Code area 8: Indexing
! ------------------------------------------------------------------------------

!!****h* Conquest/bucket_module *
!!  NAME
!! 
!!  PURPOSE
!!   this module is for definition of the derived types <local_bucket>
!!   and <remote_bucket>
!!    set maximum number
!!    allocation some variables in the above derived types
!!    deallocation some variables 
!!  
!!   For <local_bucket>, we will define in ... like 
!!   type(local_bucket) :: locbucket(mxtype) ! and mxtype=2
!!    type 1  : <phi_i|phi_j> , <phi_i|H^{local}|phi_j>
!!    type 2  : <phi_i|chi_j>
!!   On the other hand, we will use <remote_bucket>
!!   type(remote_bucket) :: rembucket(mxtype) ! and mxtype=3
!!    type 1  : <phi_i|phi_j>=Sij
!!    type 2  : <phi_i|chi_j>=Pij   
!!    type 3  : <phi_i|H|phi_j>=Hij : here H includes nonlocal parts
!! 
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!!   23/05/2000
!!  MODIFICATION HISTORY
!!   22/03/2002 dave
!!    Added RCS tags and ROBODoc header
!!   07:49, 08/01/2003 dave
!!    Added routines to deallocate bucket structures
!!   2008/02/04 08:34 dave
!!    Changed for output to file not stdout
!!   2008/05/16 ast
!!    Added timers
!!***
module bucket_module

  use datatypes
  use global_module, ONLY: io_lun
  use naba_blk_module, ONLY:naba_atm_of_blk
  use timer_stdclocks_module, ONLY: start_timer,stop_timer,tmr_std_allocation

  implicit none

  ! RCS tag for object file identification 
  character(len=80), save, private :: &
       RCSid = "$Id$"

!!****s* bucket_module/local_bucket *
!!  NAME
!!   local_bucket
!!  PURPOSE
!!   Contains information about local buckets (used when
!!   forming partial contributions to matrix elements
!!   after integration)
!!  AUTHOR
!!   T.Miyazaki
!!  SOURCE
!!
  type local_bucket
     real(double) :: rcut          ! cutoff radius
     integer      :: ind_left      ! shows index of left function
     integer      :: ind_right     ! shows index of right function
     integer      :: mx_recv_node  ! no. of remote nodes which will receive
     integer      :: no_recv_node  ! my domain partial contribution
     integer      :: no_halo_atom1 ! no. of halo atoms for a left function
     integer      :: no_halo_atom2 ! no. of halo atoms for a right function
!    integer      :: mx_pair       ! no. of pairs (mu,nu) which touch
     integer      :: no_pair       ! some of my prim blks & R(mu,nu)<rcut
     integer      :: no_pair_orb   ! some of my prim blks & R(mu,nu)<rcut
     ! the next three variables should be kept during a simulation
     integer, pointer :: list_recv_node(:) ! (mx_recv_node)
     integer, pointer :: ibegin_pair(:)    ! (mx_recv_node)
     integer, pointer :: ibegin_pair_orb(:)! (mx_recv_node)
     integer, pointer :: i_h2d(:,:)        ! (no_halo_atom1,no_halo_atom2)
  end type local_bucket
!!***

  type rem_bucket_pair
     integer :: iprim
     integer :: jhalo
     integer :: ibegin_pair
  end type rem_bucket_pair


!!****s* bucket_module/remote_bucket *
!!  NAME
!!   remote_bucket
!!  PURPOSE
!!   Contains information about remote buckets (used when 
!!   distributing partial contributions to matrix elements
!!   after integration)
!!  AUTHOR
!!   T.Miyazaki
!!  SOURCE
!!
  type remote_bucket
     type(local_bucket),pointer :: locbucket
     integer      :: mx_send_node          ! max. no. of sending nodes 
     integer      :: no_send_node          ! no. of sending nodes 
     integer      :: mx_pair_comm          ! >= max(no_of_pair(nn))
     integer,pointer :: list_send_node(:)  ! (mx_send_node)
     integer,pointer :: no_of_pair(:)      ! (mx_send_node)
     integer,pointer :: no_of_pair_orbs(:)      ! (mx_send_node)
     type(rem_bucket_pair),pointer :: bucket(:,:)        ! (mx_pair_comm,mx_send_node)
  end type remote_bucket
!!***


contains

! ------------------------------------------------------------------------------
! Subroutine set_local_bucket
! ------------------------------------------------------------------------------

!!****f* bucket_module/set_local_bucket *
!!
!!  NAME 
!!   set_local_bucket
!!  USAGE
!! 
!!  PURPOSE
!!   Sets up indexing for local buckets
!!  INPUTS
!! 
!! 
!!  USES
!!   datatypes, set_blipgrid_module
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!!   23/05/2000
!!  MODIFICATION HISTORY
!!   22/03/2002 dave
!!    Added ROBODoc headers
!!   2008/05/16 ast
!!    Added timer
!!  SOURCE
!!
!  subroutine set_local_bucket(set,rcut_in,mx_node_in,mx_pair_in,ind_left,ind_right)
  subroutine set_local_bucket(set,rcut_in,mx_node_in,ind_left,ind_right)

    use datatypes
    use set_blipgrid_module, ONLY: halo_atm
    use GenComms, ONLY: cq_abort
    use global_module, ONLY: area_index
    use memory_module, ONLY: reg_alloc_mem, type_int

    implicit none

    ! Passed variables
    type(local_bucket),intent(out):: set
    real(double),intent(in):: rcut_in
    integer,intent(in)::mx_node_in,ind_left,ind_right
    !integer,intent(in)::mx_node_in,mx_pair_in,ind_left,ind_right

    ! Local variables
    integer :: no_halo_left,no_halo_right
    integer :: stat,irc,ierr

    ! left and right functions
    set%ind_left  = ind_left
    set%ind_right = ind_right
    no_halo_left  = halo_atm(ind_left )%no_of_atom
    no_halo_right = halo_atm(ind_right)%no_of_atom

    set%rcut=rcut_in
    set%no_halo_atom1=no_halo_left ; set%no_halo_atom2=no_halo_right
    set%mx_recv_node=mx_node_in    !; set%mx_pair=mx_pair_in

    call start_timer(tmr_std_allocation)
    allocate(set%list_recv_node(set%mx_recv_node),STAT=stat)
    if(stat/=0) call cq_abort('Error alloc memory to local_bucket(list_recv_node) !')
    allocate(set%ibegin_pair(set%mx_recv_node),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to local_bucket(ibegin_pair) !')
    allocate(set%ibegin_pair_orb(set%mx_recv_node),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to local_bucket(ibegin_pair_orb) !')
    allocate(set%i_h2d(set%no_halo_atom2,set%no_halo_atom1),STAT=stat)
!    write(io_lun,*) ' allocation for loc_bucket%i_h2d', set%no_halo_atom2,set%no_halo_atom1
    if(stat/=0) call cq_abort('Error allocating memory to local_bucket(i_h2d) !')
    call reg_alloc_mem(area_index,3*set%mx_recv_node+set%no_halo_atom2*set%no_halo_atom1,type_int)
    call stop_timer(tmr_std_allocation)
    return
  end subroutine set_local_bucket
!!***

! ------------------------------------------------------------------------------
! Subroutine set_remote_bucket
! ------------------------------------------------------------------------------

!!****f* bucket_module/set_remote_bucket *
!!
!!  NAME 
!!   set_remote_bucket
!! 
!!  PURPOSE
!!   Allocates remove bucket derived type
!!  INPUTS
!! 
!! 
!!  USES
!!   GenComms
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!!   23/05/2000
!!  MODIFICATION HISTORY
!!   22/03/2002 dave
!!    Added ROBODoc headers
!!   2008/05/16 ast
!!    Added timer
!!   2009/07/08 16:43 dave
!!    Added error sizes
!!  SOURCE
!!
  subroutine set_remote_bucket(set,mx_node_in,mx_pair_in)

    use GenComms, ONLY: cq_abort
    use global_module, ONLY: area_index
    use memory_module, ONLY: reg_alloc_mem, type_int

    implicit none

    ! Passed variables
    type(remote_bucket),intent(out)::set

    integer,intent(in)::mx_node_in,mx_pair_in

    ! Local variables
    integer :: stat,irc,ierr

    set%mx_send_node=mx_node_in
    set%mx_pair_comm=mx_pair_in

    call start_timer(tmr_std_allocation)
!    write(io_lun,*) 'remote_bucket%list_send_node',set%mx_send_node
    allocate(set%list_send_node(set%mx_send_node),STAT=stat)
    if(stat/=0) call cq_abort('Err alloc memory to remote_bucket(list_send_node) ! ',set%mx_send_node)
!    write(io_lun,*) 'remote_bucket%no_of_pair',set%mx_send_node
    allocate(set%no_of_pair(set%mx_send_node),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to remote_bucket(no_of_pair) ! ',set%mx_send_node)
    allocate(set%no_of_pair_orbs(set%mx_send_node),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to remote_bucket(no_of_pair_orbs) ! ',set%mx_send_node)
!    write(io_lun,*) 'remote_bucket%bucket',set%mx_pair_comm,set%mx_send_node
    allocate(set%bucket(set%mx_pair_comm,set%mx_send_node),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to remote_bucket(bucket) !',set%mx_pair_comm,set%mx_send_node)
    call reg_alloc_mem(area_index,3*set%mx_send_node+set%mx_pair_comm*set%mx_send_node,type_int)
    call stop_timer(tmr_std_allocation)
    return
  end subroutine set_remote_bucket
!!***

! ------------------------------------------------------------------------------
! Subroutine reset_remote_bucket
! ------------------------------------------------------------------------------

!!****f* bucket_module/reset_remote_bucket *
!!
!!  NAME 
!!   reset_remote_bucket
!! 
!!  PURPOSE
!!   Reallocates remove bucket derived type
!!  INPUTS
!! 
!! 
!!  USES
!!   GenComms
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!!   23/05/2000
!!  MODIFICATION HISTORY
!!   2008/05/16 ast
!!    Added timer
!!   2009/07/08 16:43 dave
!!    Added ROBODoc headers
!!    Added error sizes
!!  SOURCE
!!
  subroutine reset_remote_bucket(set)

    use GenComms, ONLY: cq_abort

    implicit none

    ! Passed variables
    type(remote_bucket),intent(inout)::set

    integer :: stat

    call start_timer(tmr_std_allocation)
    deallocate(set%bucket,STAT=stat)
    allocate(set%bucket(set%mx_pair_comm,set%mx_send_node),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to remote_bucket(bucket) !',set%mx_pair_comm,set%mx_send_node)
    call stop_timer(tmr_std_allocation)
    return
  end subroutine reset_remote_bucket
!!***

! ------------------------------------------------------------------------------
! Subroutine free_local_bucket
! ------------------------------------------------------------------------------

!!****f* bucket_module/free_local_bucket *
!!
!!  NAME 
!!   free_local_bucket
!!  USAGE
!! 
!!  PURPOSE
!!   Frees up local bucket memory
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   2003/01/08
!!  MODIFICATION HISTORY
!!   2008/05/16 ast
!!    Added timer
!!  SOURCE
!!
  subroutine free_local_bucket(set)

    use GenComms, ONLY: cq_abort
    use global_module, ONLY: area_index
    use memory_module, ONLY: reg_dealloc_mem, type_int

    implicit none

    ! Passed variables
    type(local_bucket), intent(inout) :: set

    ! Local variables
    integer :: stat

    call reg_dealloc_mem(area_index,3*set%mx_recv_node+set%no_halo_atom2*set%no_halo_atom1,type_int)
    ! Zero values
    set%ind_left  = 0
    set%ind_right = 0
    set%rcut=0.0_double
    set%no_halo_atom1=0 ; set%no_halo_atom2=0
    set%mx_recv_node=0    !; set%mx_pair=0
    ! Deallocate memory
    call start_timer(tmr_std_allocation)
    deallocate(set%i_h2d,set%ibegin_pair_orb,set%ibegin_pair,set%list_recv_node,STAT=stat)
    if(stat/=0) call cq_abort('Error deallocating memory from local_bucket !')
    call stop_timer(tmr_std_allocation)
    return
  end subroutine free_local_bucket
!!***

! ------------------------------------------------------------------------------
! Subroutine free_remote_bucket
! ------------------------------------------------------------------------------

!!****f* bucket_module/free_remote_bucket *
!!
!!  NAME 
!!   free_remote_bucket
!!  USAGE
!! 
!!  PURPOSE
!!   Frees up remote bucket memory
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   2003/01/08
!!  MODIFICATION HISTORY
!!   2008/05/16 ast
!!    Added timer
!!  SOURCE
!!
  subroutine free_remote_bucket(set)

    use GenComms, ONLY: cq_abort
    use memory_module, ONLY: reg_dealloc_mem, type_int
    use global_module, ONLY: area_index

    implicit none

    ! Passed variables
    type(remote_bucket), intent(inout) :: set

    ! Local variables
    integer :: stat

    call start_timer(tmr_std_allocation)
    call reg_dealloc_mem(area_index,3*set%mx_send_node+set%mx_pair_comm*set%mx_send_node,type_int)
    ! Deallocate memory
    deallocate(set%bucket,set%no_of_pair_orbs,set%no_of_pair,set%list_send_node,STAT=stat)
    if(stat/=0) call cq_abort('Error deallocating memory from remote_bucket !')
    call stop_timer(tmr_std_allocation)
    return
  end subroutine free_remote_bucket
!!***

end module bucket_module
