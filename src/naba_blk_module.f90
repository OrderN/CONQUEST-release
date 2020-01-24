! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module naba_blk_module
! ------------------------------------------------------------------------------
! Code area 12: integration
! ------------------------------------------------------------------------------

!!****h* Conquest/naba_blk_module
!!  NAME
!!   naba_blk_module
!!  PURPOSE
!!   this module includes the definition of the following derived types
!!     1. naba_blk_of_atm  : neighbour blocks of primary atoms
!!     2. naba_atm_of_blk  : neighbour atoms of primary blocks
!!     3. halo_atm_of_blk  : halo atoms of primary blocks
!!     4. comm_in_BtoG                                    
!! 
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!!   31/05/2000
!!  MODIFICATION HISTORY
!!   09:50, 29/10/2002 dave
!!    Added ROBODoc headers
!!   13:52, 10/02/2003 drb 
!!    Added deallocation routines (or at least tidied them)
!!   2008/05/16 ast
!!    Added timers
!!   2014/09/15 18:30 lat
!!    fixed call start/stop_timer to timer_module (not timer_stdlocks_module !)
!!   2016/07/06 17:30 nakata
!!    Renamed comm_in_BG -> comm_in_BtoG
!!  SOURCE
!!
module naba_blk_module

  use datatypes
  use timer_module,           only: start_timer,stop_timer
  use timer_stdclocks_module, only: tmr_std_allocation

  implicit none

  !---------------------------------------------------------------------
  !1. naba blocks for each primary atom
  !---------------------------------------------------------------------
  type naba_blk_of_atm
     ! mx_iprim         : = max{parts%mx_iprim} within all nodes
     ! mx_naba_blk      : >= max{no_naba_blk(iprim)} 
     ! no_naba_blk      : no. of naba blks for each primary atom (iprim)
     ! nxmin etc.       : min of x-coord of naba blks for each iprim
     !                     will be used in BlipGrid transforms
     ! list_naba_blk    : CC label of neighbour blocks in cov. sets
     ! send_naba_blk    : CC label of neighbour blocks in sim. cell
     ! offset_naba_blk  : offset of neighbour blocks
     !                     calculated in ncover's of remote DCS_parts
     integer :: mx_iprim,mx_naba_blk
     integer,pointer:: no_naba_blk(:)  ! (mx_iprim)
     integer,pointer:: nxmin(:)        ! (mx_iprim)
     integer,pointer:: nxmax(:)        ! (mx_iprim)
     integer,pointer:: nymin(:)        ! (mx_iprim)
     integer,pointer:: nymax(:)        ! (mx_iprim)
     integer,pointer:: nzmin(:)        ! (mx_iprim)
     integer,pointer:: nzmax(:)        ! (mx_iprim)
     integer,pointer:: list_naba_blk(:,:)   ! (mx_naba_blk,mx_iprim)
     integer,pointer:: send_naba_blk(:,:)   ! (mx_naba_blk,mx_iprim)
     integer,pointer:: offset_naba_blk(:,:) ! (mx_naba_blk,mx_iprim)
  end type naba_blk_of_atm

  !---------------------------------------------------------------------
  !2. naba atoms for each primary block
  !---------------------------------------------------------------------
  type naba_atm_of_blk
     integer :: mx_iprim_blk  ! max. no. of primary blocks = domain%mx_ngonn
     integer :: mx_part       ! max. no. of naba parts for each block
     integer :: mx_atom       ! max. no. of naba atoms for each block
     integer,pointer :: no_of_part(:)          ! (mx_iprim_blk)
     integer,pointer :: no_of_atom(:)          ! (mx_iprim_blk)
     integer,pointer :: ibegin_blk(:)          ! (mx_iprim_blk)
     integer,pointer :: list_part(:,:)         ! (mx_part,mx_iprim_blk)
     integer,pointer :: no_atom_on_part(:,:)   ! (mx_part,mx_iprim_blk)
     integer,pointer :: ibegin_part(:,:)       ! (mx_part,mx_iprim_blk)
     integer,pointer :: list_atom(:,:)         ! (mx_atom,mx_iprim_blk)
     integer,pointer :: list_atom_by_halo(:,:) ! (mx_atom,mx_iprim_blk)
!TM VARNSF : START
     integer :: function_type  ! shows the type of functions
     integer, pointer :: no_of_orb(:)      ! (mx_iprim_blk)  number of naba orbitals for each primary block
     integer, pointer :: ibegin_blk_orb(:)     ! (mx_iprim_blk)
     integer, pointer :: ibeg_orb_atom(:,:)  ! (mx_atom,mx_iprim_blk) 
!TM VARNSF : END
     !no_of_part       : no. of naba parts for each primary block
     !no_of_atom       : no. of naba atoms for each primary block
     !ibegin_blk       : seq. no. of 1st naba atoms for each primary block
     !list_part        : labels of naba parts (NOPG in DCS_parts)
     !no_atom_on_part  : no. of naba atoms in each naba part
     !ibegin_part      : seq. no. of the first naba atm for each prim. atom
     !list_atom        : parts seq. labels of naba atms
     !list_atom_by_halo: halo seq. labels of naba atms
  end type naba_atm_of_blk

  !---------------------------------------------------------------------
  !3. halo atoms and partitions of my domain
  !---------------------------------------------------------------------
  type halo_atm_of_blk
     !type(naba_atm_of_blk),pointer :: naba_atoms_of_blocks
     integer :: mx_node
     integer :: mx_part         !
     integer :: mx_mcover       ! = DCS_parts%mx_mcover
     integer :: no_of_node
     integer :: no_of_part
     integer :: no_of_atom
     integer,pointer :: list_of_node(:) ! (mx_node)
     integer,pointer :: list_of_part(:) ! (mx_part)
     integer,pointer :: ihalo(:)        ! (mx_mcover)
!TM VARNSF : START
     integer, pointer :: norb(:)        ! (no_of_atom)   ?? OR global id, OR species ??
!TM VARNSF : END
     !no_of_node       : no. of nodes responsible for the halo atoms
     !no_of_part       : no. of halo parts which includes halo atoms
     !no_of_atom       : no. of halo atoms
     !list_of_node(:)  : this will be used only for checking
     !list_of_part(:)  : NOPG label of halo partitions
     !ihalo(:)         : transcription table (cover -> halo atoms)
     !  We need this to know the index of a pair of two atoms (cover seq.
     !  label) through local_bucket%i_h2d(halo1,halo2) in calculating
     !  domain partial contributions of matrix elements.
  end type halo_atm_of_blk

  !---------------------------------------------------------------------
  !4. Variables used in the communication on BlipGrid transformations
  !---------------------------------------------------------------------
  type comm_in_BtoG
     integer :: mx_iprim    !   = bundle%mx_iprim
     integer :: mx_recv_node,mx_send_node
     integer :: mx_pair     ! >= max[no_sent_pairs(:,:)]
     integer :: mx_recv_call! >= max[no_recv_call(:,:)]
     !for Bundle responsible node (for sending)
     integer,pointer :: no_recv_node(:)      !(mx_iprim)
     integer,pointer :: list_recv_node(:,:)  !(mx_recv_node,mx_iprim)
     ! integer,pointer :: ibegin_naba_blk(:,:) !(mx_recv_node,mx_iprim) **********
     integer,pointer :: no_naba_blk(:,:)     !(mx_recv_node,mx_iprim) **********
     !for Domain responsible node (for receiving)
     integer,pointer :: no_send_node(:)      !(mx_iprim)
     integer,pointer :: list_send_node(:,:)  !(mx_send_node,mx_iprim)
     integer,pointer :: no_sent_pairs(:,:)   !(mx_send_node,mx_iprim)
     !table showing where to put sent BtoG transformed functions
     integer,pointer :: table_blk(:,:)       !(mx_pair,mx_recv_call)
     integer,pointer :: table_pair(:,:)      !(mx_pair,mx_recv_call)
     integer,pointer :: ibeg_recv_call(:)      !(mx_iprim)
     !Blip-grid transformations are done one by one with respect to
     !primary atoms. Thus, most of the above arrays use the final 
     !rank(?) to identify a primary number of atoms.
     !
     !For tables, a serial no. of receiving calls is used to identify
     !pairs of an atom and its neighbour blocks to reduce the memory
     !size. (The arrays like table_blk(mx_pair,mx_send_node,mx_iprim) 
     !might be easier to understand, but needs much more memory.)
     !mx_iprim     : max. no. of primary atoms for all bundles
     !no_recv_node : no. of receiving nodes in sending each BtoG 
     !               transformed support function
     !mx_recv_node : maximum of no_recv_node
     !no_send_node : no. of sending nodes in receiving each BtoG 
     !               transformed support function
     !mx_send_node : maximum of no_send_node
     !no_sent_pairs: no. of pairs(an atom and its naba blks) for 
     !               each receiving call
     !mx_pair      : maximum of no_sent_pairs
     !(no_recv_call : no. of receiving calls for each primary seq. no.
     ! this has been replaced by ibeg_recv_call)
     !ibeg_recv_call: First seq. no. of receiving calls for each primary atom no.
     !mx_recv_call : maximum of no_recv_call
     !table_blk    : shows the no. of a primary block corresponding 
     !               to each pair sent from a Bundle responsible node
     !               to my node(Domain Responsible node).
     !table_pair   : shows a index of a naba atom corresponding
     !               to each pair sent from ... 
     !               (The neighbour atom is identified through the naba
     !               list for one of my primary blocks)
  end type comm_in_BtoG

  !!***

contains

!!****f* naba_blk_module/alloc_naba_blk *
!!
!!  NAME 
!!   alloc_naba_blk
!!  USAGE
!! 
!!  PURPOSE
!! 
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!!   July 2000
!!  MODIFICATION HISTORY
!!   2008/05/16 ast
!!    Added timers
!!  SOURCE
!!
  subroutine alloc_naba_blk(set,mx1,mx2)

    use GenComms, ONLY: cq_abort
    use global_module, ONLY: area_integn
    use memory_module, ONLY: reg_alloc_mem, type_int

    implicit none

    integer,intent(in)::mx1,mx2
    type(naba_blk_of_atm),intent(out)::set
    integer :: stat

    call start_timer(tmr_std_allocation)
    set%mx_iprim=mx1 ; set%mx_naba_blk=mx2
    allocate(set%no_naba_blk(mx1),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to no_naba_blk !',stat)
    allocate(set%nxmin(mx1),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to nxmin !',stat)
    allocate(set%nxmax(mx1),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to nxmax !',stat)
    allocate(set%nymin(mx1),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to nymin !',stat)
    allocate(set%nymax(mx1),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to nymax !',stat)
    allocate(set%nzmin(mx1),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to nzmin !',stat)
    allocate(set%nzmax(mx1),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to nzmax !',stat)
    allocate(set%list_naba_blk(mx2,mx1),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to list_naba_blk !',stat)
    allocate(set%send_naba_blk(mx2,mx1),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to send_naba_blk !',stat)
    allocate(set%offset_naba_blk(mx2,mx1),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to offset_naba_blk !',stat)
    call reg_alloc_mem(area_integn, 7*mx1+3*mx1*mx2,type_int)
    call stop_timer(tmr_std_allocation)
    return
  end subroutine alloc_naba_blk
!!***

!!****f* naba_blk_module/dealloc_naba_blk *
!!
!!  NAME 
!!   dealloc_naba_blk
!!  USAGE
!! 
!!  PURPOSE
!! 
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!!   July 2000
!!  MODIFICATION HISTORY
!!   2006/10/20 09:56 dave
!!    Changed attribute of set
!!   2008/05/16 ast
!!    Added timers
!!  SOURCE
!!
  subroutine dealloc_naba_blk(set)

    use GenComms, ONLY: cq_abort
    use global_module, ONLY: area_integn
    use memory_module, ONLY: reg_dealloc_mem, type_int

    implicit none

    type(naba_blk_of_atm),intent(inout)::set
    integer :: stat

    call start_timer(tmr_std_allocation)
    call reg_dealloc_mem(area_integn, 7*set%mx_iprim+3*set%mx_iprim*set%mx_naba_blk,type_int)
    deallocate(set%offset_naba_blk, set%send_naba_blk, set%list_naba_blk,&
         set%nzmax, set%nzmin, set%nymax, set%nymin, set%nxmax, set%nxmin,&
         set%no_naba_blk, STAT=stat)
    if(stat /= 0) call cq_abort(' ERROR in deallocating naba_blk ',stat)
    call stop_timer(tmr_std_allocation)
    return
  end subroutine dealloc_naba_blk
!!***

!!****f* naba_blk_module/alloc_naba_atm *
!!
!!  NAME 
!!   alloc_naba_atm
!!  USAGE
!! 
!!  PURPOSE
!! 
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!!   July 2000
!!  MODIFICATION HISTORY
!!   2006/07/04 08:14 dave
!!    Added allocates for variable NSF arrays
!!   2008/05/16 ast
!!    Added timers
!!   2016/07/20 16:30 nakata
!!    Comment: alloc_naba_atm is to allocate naba_atoms_of_blocks
!!             (renamed from naba_atm, passed in "set")
!!  SOURCE
!!
  subroutine alloc_naba_atm(set,mx1,mx2,mx3)

    use GenComms, ONLY: cq_abort
    use global_module, ONLY: area_integn
    use memory_module, ONLY: reg_alloc_mem, type_int

    implicit none

    integer,intent(in)::mx1,mx2,mx3
    type(naba_atm_of_blk),intent(out)::set
    integer :: stat

    call start_timer(tmr_std_allocation)
    set%mx_iprim_blk=mx1 
    set%mx_part=mx2
    set%mx_atom=mx3
    allocate(set%no_of_part(mx1),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to no_of_part !',stat)
    allocate(set%no_of_atom(mx1),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to no_of_atom !',stat)
    allocate(set%no_of_orb(mx1),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to no_of_orb !',stat)
    allocate(set%ibegin_blk(mx1),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to ibegin_blk !',stat)
    allocate(set%ibegin_blk_orb(mx1),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to ibegin_blk_orb !',stat)
    allocate(set%list_part(mx2,mx1),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to list_part !',stat)
    allocate(set%no_atom_on_part(mx2,mx1),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to no_atom_on_part !',stat)
    allocate(set%ibegin_part(mx2,mx1),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to ibegin_part !',stat)
    allocate(set%ibeg_orb_atom(mx3,mx1),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to ibeg_orb_atom !',stat)
    allocate(set%list_atom(mx3,mx1),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to list_atom !',stat)
    allocate(set%list_atom_by_halo(mx3,mx1),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to list_atom_by_halo !',stat)
    call reg_alloc_mem(area_integn,5*mx1+3*mx2*mx1+3*mx3*mx1,type_int)
    call stop_timer(tmr_std_allocation)
    return
  end subroutine alloc_naba_atm
!!***

!!****f* naba_blk_module/dealloc_naba_atm *
!!
!!  NAME 
!!   dealloc_naba_atm
!!  USAGE
!! 
!!  PURPOSE
!! 
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!!   July 2000
!!  MODIFICATION HISTORY
!!   2006/07/04 08:13 dave
!!    Added deallocate for new variable NSF arrays
!!   2006/10/20 09:56 dave
!!    Changed attribute of set
!!   2008/05/16 ast
!!    Added timers
!!   2016/07/20 16:30 nakata
!!    Comment: alloc_naba_atm is to allocate naba_atoms_of_blocks
!!             (renamed from naba_atm, passed in "set")
!!  SOURCE
!!
  subroutine dealloc_naba_atm(set)

    use GenComms, ONLY: cq_abort
    use global_module, ONLY: area_integn
    use memory_module, ONLY: reg_dealloc_mem, type_int

    implicit none

    type(naba_atm_of_blk),intent(inout)::set
    integer :: stat

    call start_timer(tmr_std_allocation)
    call reg_dealloc_mem(area_integn,5*set%mx_iprim_blk+3*set%mx_part*set%mx_iprim_blk+&
         3*set%mx_atom*set%mx_iprim_blk,type_int)
    deallocate(set%list_atom_by_halo, set%list_atom, set%ibeg_orb_atom, set%ibegin_part, &
         set%no_atom_on_part, set%list_part, set%ibegin_blk, set%ibegin_blk_orb, &
         set%no_of_orb, set%no_of_atom, set%no_of_part, STAT=stat)
    if(stat /= 0) call cq_abort('ERROR in deallocating naba_atm',stat)
    call stop_timer(tmr_std_allocation)
    return
  end subroutine dealloc_naba_atm
!!***

!!****f* naba_blk_module/alloc_halo_atm *
!!
!!  NAME 
!!   alloc_halo_atm
!!  USAGE
!! 
!!  PURPOSE
!! 
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!!   July 2000
!!  MODIFICATION HISTORY
!!   2008/05/16 ast
!!    Added timers
!!   2016/07/20 16:30 nakata
!!    Comment: alloc_halo_atm is to allocate halo_atoms_of_blocks
!!             (renamed from halo_atm, passed in "set")
!!  SOURCE
!!
  subroutine alloc_halo_atm(set,mx1,mx2,mx3)

    use GenComms, ONLY: cq_abort
    use global_module, ONLY: area_integn
    use memory_module, ONLY: reg_alloc_mem, type_int

    implicit none

    integer,intent(in)::mx1,mx2,mx3
    type(halo_atm_of_blk),intent(out)::set
    integer :: stat

    call start_timer(tmr_std_allocation)
    set%mx_node=mx1 
    set%mx_part=mx2
    set%mx_mcover=mx3

    allocate(set%list_of_node(mx1),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to list_of_node !',stat)
    allocate(set%list_of_part(mx2),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to list_of_part !',stat)
    allocate(set%ihalo(mx3),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to ihalo !',stat)
    allocate(set%norb(mx3),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to norb !',stat)
    call reg_alloc_mem(area_integn,mx1+mx2+2*mx3,type_int)
    call stop_timer(tmr_std_allocation)
    return
  end subroutine alloc_halo_atm
!!***

!!****f* naba_blk_module/dealloc_halo_atm *
!!
!!  NAME 
!!   dealloc_halo_atm
!!  USAGE
!! 
!!  PURPOSE
!! 
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!!   July 2000
!!  MODIFICATION HISTORY
!!   2006/10/20 09:56 dave
!!    Changed attribute of set
!!   2008/05/16 ast
!!    Added timers
!!   2016/07/20 16:30 nakata
!!    Comment: dealloc_halo_atm is to deallocate halo_atoms_of_blocks
!!             (renamed from halo_atm, passed in "set")
!!  SOURCE
!!
  subroutine dealloc_halo_atm(set)

    use GenComms, ONLY: cq_abort
    use global_module, ONLY: area_integn
    use memory_module, ONLY: reg_dealloc_mem, type_int

    implicit none

    type(halo_atm_of_blk),intent(inout)::set
    integer :: stat

    call start_timer(tmr_std_allocation)
    call reg_dealloc_mem(area_integn,set%mx_node+set%mx_part+2*set%mx_mcover,type_int)
    deallocate(set%norb,set%ihalo, set%list_of_part, set%list_of_node, STAT=stat)
    if(stat /= 0) call cq_abort('ERROR in deallocating halo_atm',stat)
    call stop_timer(tmr_std_allocation)
    return
  end subroutine dealloc_halo_atm
!!***

!!****f* naba_blk_module/alloc_comm_in_BtoG1 *
!!
!!  NAME 
!!   alloc_comm_in_BtoG1
!!  USAGE
!! 
!!  PURPOSE
!! 
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   T.Miyazaki/D.R.Bowler
!!  CREATION DATE
!!   July 2000 and 2006/09/05
!!  MODIFICATION HISTORY
!!   2008/05/16 ast
!!    Added timers
!!   2016/07/06 17:30 nakata
!!    Renamed subroutine alloc_comm_in_BG1 -> alloc_comm_in_BtoG1
!!    Renamed comm_in_BG -> comm_in_BtoG
!!  SOURCE
!!
  subroutine alloc_comm_in_BtoG1(set,mx1,mx2)

    use GenComms, ONLY: cq_abort

    implicit none

    integer,intent(in)::mx1,mx2
    type(comm_in_BtoG),intent(inout) :: set
    integer :: stat

    set%mx_iprim=mx1
    set%mx_recv_node=mx2

    call start_timer(tmr_std_allocation)
    allocate(set%no_recv_node(mx1),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to no_recv_node !',stat)
    allocate(set%list_recv_node(mx2,mx1),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to list_recv_node !',mx2,mx1)
    allocate(set%no_naba_blk(mx2,mx1),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to no_naba_blk !',stat)
    allocate(set%no_send_node(mx1),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to no_send_node !',stat)
    allocate(set%ibeg_recv_call(mx1),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to ibeg_recv_call !',stat)
    call stop_timer(tmr_std_allocation)
    return
  end subroutine alloc_comm_in_BtoG1
!!***

!!****f* naba_blk_module/alloc_comm_in_BtoG2 *
!!
!!  NAME 
!!   alloc_comm_in_BtoG2
!!  USAGE
!! 
!!  PURPOSE
!! 
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   T.Miyazaki/D.R.Bowler
!!  CREATION DATE
!!   July 2000 and 2006/09/05
!!  MODIFICATION HISTORY
!!   2008/05/16 ast
!!    Added timers
!!   2012/06/17 L.Tong
!!   - Changed set from intent(out) to intent(inout). Otherwise
!!     member mx_iprim becomes undefined upon entering the subroutine,
!!     and is not defined anywhere within the subroutine.
!!   2016/07/06 17:30 nakata
!!    Renamed subroutine alloc_comm_in_BG2 -> alloc_comm_in_BtoG2
!!    Renamed comm_in_BG -> comm_in_BtoG
!!  SOURCE
!!
  subroutine alloc_comm_in_BtoG2(set,mx3,mx4,mx5)

    use GenComms, ONLY: cq_abort

    implicit none

    integer,intent(in)::mx3,mx4,mx5
    type(comm_in_BtoG),intent(inout) :: set
    integer :: stat

    set%mx_send_node=mx3
    set%mx_pair     =mx4
    set%mx_recv_call=mx5

    call start_timer(tmr_std_allocation)
    allocate(set%list_send_node(mx3,set%mx_iprim),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to list_send_node !',stat)
    allocate(set%no_sent_pairs(mx3,set%mx_iprim),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to no_sent_pairs !',stat)
    allocate(set%table_blk(mx4,mx5),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to table_blk !',stat)
    allocate(set%table_pair(mx4,mx5),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to table_pair !',stat)
    call stop_timer(tmr_std_allocation)
    return
  end subroutine alloc_comm_in_BtoG2
!!***

!!****f* naba_blk_module/dealloc_comm_in_BtoG *
!!
!!  NAME 
!!   dealloc_comm_in_BtoG
!!  USAGE
!! 
!!  PURPOSE
!! 
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!!   July 2000
!!  MODIFICATION HISTORY
!!   2006/10/20 09:56 dave
!!    Changed attribute of set
!!   2008/05/16 ast
!!    Added timers
!!   2016/07/06 17:30 nakata
!!    Rename subroutine dealloc_comm_in_BG -> dealloc_comm_in_BtoG
!!    Renamed comm_in_BG -> comm_in_BtoG and comBG to comm_naba_blocks_of_atoms
!!  SOURCE
!!
  subroutine dealloc_comm_in_BtoG(set)

    use GenComms, ONLY: cq_abort

    implicit none

    type(comm_in_BtoG),intent(inout) :: set
    integer :: stat

    call start_timer(tmr_std_allocation)
    deallocate(set%ibeg_recv_call, set%table_pair, set%table_blk, &
         set%no_sent_pairs, set%list_send_node, set% no_send_node, &
         set%no_naba_blk, set%list_recv_node, &
         set%no_recv_node, STAT=stat)
    if(stat /= 0) call cq_abort(' ERROR in deallocating comm_naba_blocks_of_atoms ',stat)
    call stop_timer(tmr_std_allocation)
    return
  end subroutine dealloc_comm_in_BtoG
!!***

end module naba_blk_module
