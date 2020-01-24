! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module set_bucket_module
! ------------------------------------------------------------------------------
! Code area 8: indexing
! ------------------------------------------------------------------------------

!!****h* Conquest/set_bucket_module *
!!  NAME
!! 
!!  PURPOSE
!!     This module includes the arrays (local and remote buckets) 
!!    and subroutines which makes these arrays.
!!
!!     Local bucket will be used in summing up block partial contributions
!!    to the local array (domain partial contribution) of the matrix 
!!    elements and send them.
!!
!!     Remote Bucket provides information of sending nodes, where to put
!!    each sent domain partial contributions into the local array of 
!!    the matrix
!!
!!   Structure of the subroutines in this module 
!!    set_bucket
!!     |--  (set_remote_bucket),(set_local_bucket)  <- bucket_module
!!   (contains)
!!     |--   make_sendinfo_ME  making the information of sending nodes
!!     |--   make_pair_DCSpart makes pairs of atoms, which has a non-zero
!!     |                       domain partial contribution of my node.
!!     |                       also makes the information of recv nodes
!!     |--   setup_sendME      sends the information of pairs made in 
!!     |                       make_pair_DCSpart
!!     |--   recv_npairME      receives rem_bucket%no_of_pair
!!     |--   setup_recvME      receives the infomation of pairs.
!!                               makes rem_bucket%bucket, which shows how to
!!                             put the received arrays into the local arrays
!! 
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!!   30/05/2000
!!  MODIFICATION HISTORY
!!   11/05/2001 drb
!!    Changed bundle%iprim_group to BCS_parts%iprim_group - is it the right CS ? 
!!   14:51, 26/02/2003 drb 
!!    Added header and tweaked a use
!!   11:41, 10/07/2006 drb 
!!    Added new buckets for PAO integrals on grid
!!   2006/10/18 17:20 dave
!!    Changes from TM incorporated to find max recvs/sends
!!   2008/02/06 08:36 dave
!!    Changed for output to file not stdout
!!   2008/05/16
!!    Added timers
!!   2014/09/15 18:30 lat
!!    fixed call start/stop_timer to timer_module (not timer_stdlocks_module !)
!!   2016/07/15 18:30 nakata
!!    Renamed sf_sf_loc -> atomf_atomf_loc, sf_nlpf_loc -> atomf_nlpf_loc, 
!!            sf_sf_rem -> atomf_atomf_rem, sf_nlpf_rem -> atomf_nlpf_rem, sf_H_sf_rem -> atomf_H_atomf_rem
!!   2016/12/29 19:00 nakata
!!    Removed pao_sf_loc and pao_H_sf_rem, which are no longer needed 
!!    Instead, we will use atomf_atomf_loc*C(sf,atomf) and atomf_H_atomf_rem*C(sf,atomf)
!!  SOURCE
!!
module set_bucket_module

  use bucket_module
  use global_module,          only: io_lun
  use numbers,                only: RD_ERR
  use matrix_module,          only: matrix_halo
  use timer_module,           only: start_timer,stop_timer
  use timer_stdclocks_module, only: tmr_std_indexing

  implicit none

  integer,parameter  :: mx_type_loc=3
  integer,parameter  :: mx_type_rem=4
  type(local_bucket),target :: loc_bucket(mx_type_loc)
  type(remote_bucket) :: rem_bucket(mx_type_rem)

  ! Parameters for indexing buckets
  integer, parameter :: atomf_atomf_loc = 1
  integer, parameter :: atomf_nlpf_loc  = 2
  integer, parameter :: pao_sf_loc      = 3
  integer, parameter :: atomf_atomf_rem    = 1
  integer, parameter :: atomf_nlpf_rem     = 2
  integer, parameter :: atomf_H_atomf_rem  = 3

!!***

contains

!!****f* set_bucket_module/set_bucket *
!!
!!  NAME 
!!   set_bucket
!!  USAGE
!! 
!!  PURPOSE
!!   main subroutine in this module
!!   contains other subroutines to share the arrays.
!! (isend_array,recv_array,nsend_req,nrecv_stat etc.)
!!
!!   makes local_bucket and remote_bucket
!!
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!!   May 2000
!!  MODIFICATION HISTORY
!!   2016/07/06 17:30 nakata
!!    Renamed comBG -> comm_naba_blocks_of_atoms
!!   2016/07/15 18:30 nakata
!!    Renamed sf_sf_loc -> atomf_atomf_loc, sf_nlpf_loc -> atomf_nlpf_loc,
!!            sf_sf_rem -> atomf_atomf_rem, sf_nlpf_rem -> atomf_nlpf_rem, sf_H_sf_rem -> atomf_H_atomf_rem
!!   2016/08/01 17:30 nakata
!!    Introduced atomf instead of sf and paof
!!   2016/12/29 19:30 nakata
!!    Removed max_PAOPpair, which is no longer needed 
!!
!!  SOURCE
!!
  subroutine set_bucket(myid)

    use datatypes
    use global_module, ONLY: iprint_index, sf, nlpf, paof, atomf, numprocs, &
                             flag_basis_set, PAOs, blips
    use bucket_module
!    use maxima_module, ONLY: mx_send_node_BtoG, mx_recv_node_BtoG, &
!         mx_Spair_remnode,  mx_Spair_domain,  &
!         mx_APpair_remnode, mx_APpair_domain
    use set_blipgrid_module,ONLY: comm_naba_blocks_of_atoms
    use comm_array_module,  ONLY: isend_array,irecv_array!,check_commarray_int!,check_commarray_real
    use matrix_data, ONLY: halo, rcut, aSa_range, APrange, aHa_range
    use GenComms, ONLY: my_barrier, cq_abort
    use mpi, ONLY: MPI_STATUS_SIZE

    implicit none

    integer :: nsend_req(3*numprocs)
    integer :: nrecv_stat(MPI_STATUS_SIZE,3*numprocs)
    integer,intent(in) :: myid
    integer :: myid_ibegin,myid_npair, myid_npair_orb
    integer :: ii, stat
    integer :: inode,mx_pair1,mx_pair2,mx_pair3,mx_pair_orb1,mx_pair_orb2,mx_pair_orb3
    real(double) :: rcut_S,rcut_SP
    integer :: max_Spair, max_APpair
    integer :: mx_recv_BRnode ! for local bucket
    integer :: mx_send_DRnode ! for remote bucket
    !Automatic Array
    integer :: buff_npair(numprocs), buff_npair_orb(numprocs)

    call start_timer(tmr_std_indexing)

    !cutoff radii for matrices
    rcut_S = rcut(aSa_range)
    rcut_SP= rcut(APrange)

    !atomf = sf (for blips and one_to_one PAOs) or paof (for contracted PAOs)
    !For type 1 of local_bucket   : Sij=<atomf_i|atomf_j>, and H^(local)_ij
    !For type 2 of local_bucket   : Pij=<atomf_i|chi_j>,|chi_j> proj. func.

    !For type 1 of remote_bucket  : Sij=<atomf_i|atomf_j>
    !For type 2 of remote_bucket  : Pij=<atomf_i|chi_j>,|chi_j> proj. func.
    !For type 3 of remote_bucket  : Hij=<atomf_i|H|atomf_j>
    ! Here H includes nonlocal parts

    !calculate mx_send_node
    call calc_mx_nodes_bucket(myid, comm_naba_blocks_of_atoms, mx_recv_BRnode, mx_send_DRnode)

    call set_local_bucket(loc_bucket(atomf_atomf_loc),rcut_S, mx_recv_BRnode,atomf,atomf)
    call set_local_bucket(loc_bucket(atomf_nlpf_loc),rcut_SP, mx_recv_BRnode,atomf,nlpf)
    !** type 1 local_bucket -> type 1&3 remote_bucket **
    call my_barrier()
    call make_pair_DCSpart(loc_bucket(atomf_atomf_loc),max_Spair)

    ! allocate & set  remote bucket
    !   Now we assume the left functions are always support functions.
    !  Thus, no. of domain responsible nodes, which send their domain 
    !  partial contributions of matrix elements, is no. of receiving nodes
    !  in the communication of Blip-Grid transformations.
    ! 2006 Oct 18 by TM
    !  Now, the number of nodes for local and remote buckets are calculated
    !  in subroutine 'calc_mx_nodes_bucket'.

    rem_bucket(atomf_atomf_rem)%locbucket => loc_bucket(atomf_atomf_loc)
    rem_bucket(atomf_nlpf_rem)%locbucket => loc_bucket(atomf_nlpf_loc)
    rem_bucket(atomf_H_atomf_rem)%locbucket => loc_bucket(atomf_atomf_loc)

    call set_remote_bucket(rem_bucket(atomf_atomf_rem),mx_send_DRnode,max_Spair)
    call set_remote_bucket(rem_bucket(atomf_H_atomf_rem),mx_send_DRnode,max_Spair)

    ! Bundle responsible node makes rem_bucket%no_send_node and 
    !list_send_node to prepare receiving the information from its
    !remote nodes (domain responsible nodes).
    ! WARNING!!
    ! Now, we assume the left functions are always support functions.
    call make_sendinfo_ME(rem_bucket(atomf_atomf_rem),comm_naba_blocks_of_atoms)
    rem_bucket(atomf_H_atomf_rem)%no_send_node=rem_bucket(atomf_atomf_rem)%no_send_node
    rem_bucket(atomf_H_atomf_rem)%list_send_node(:)=rem_bucket(atomf_atomf_rem)%list_send_node(:)
    !Domain responsible node sends its information of pairs of naba atms
    call setup_sendME(myid,myid_ibegin,myid_npair,myid_npair_orb,loc_bucket(atomf_atomf_loc))
    !Bundle responsible node receives rem_bucket%no_of_pair and no_of_pair_orb
    call recv_npairME(myid,myid_npair,myid_npair_orb,rem_bucket(atomf_atomf_rem))
    rem_bucket(atomf_H_atomf_rem)%no_of_pair(:) = rem_bucket(atomf_atomf_rem)%no_of_pair(:)
    rem_bucket(atomf_H_atomf_rem)%no_of_pair_orbs(:) = rem_bucket(atomf_atomf_rem)%no_of_pair_orbs(:)
    if(rem_bucket(atomf_H_atomf_rem)%mx_pair_comm<rem_bucket(atomf_atomf_rem)%mx_pair_comm) then
       rem_bucket(atomf_H_atomf_rem)%mx_pair_comm=rem_bucket(atomf_atomf_rem)%mx_pair_comm
       call reset_remote_bucket(rem_bucket(atomf_H_atomf_rem))
    end if
    !Bundle responsible node receives the information of pairs of naba
    ! atoms of its remote nodes and Constructs RemoteBucket
    call setup_recvME(atomf_atomf_rem,myid,myid_ibegin,rem_bucket(atomf_atomf_rem),halo(aSa_range), &
         rem_bucket(atomf_H_atomf_rem),halo(aHa_range))
    deallocate(isend_array,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating isend_array in set_bucket: ",stat)

    !** type 2 local_bucket -> type 2 remote_bucket **
    ! same as above (dummy arg. for make_pair_DCSpart is different)

    call make_pair_DCSpart(loc_bucket(atomf_nlpf_loc),max_APpair)
    call set_remote_bucket(rem_bucket(atomf_nlpf_rem),mx_send_DRnode,max_APpair)
    rem_bucket(atomf_nlpf_rem)%no_send_node=rem_bucket(atomf_atomf_rem)%no_send_node
    rem_bucket(atomf_nlpf_rem)%list_send_node(:)=rem_bucket(atomf_atomf_rem)%list_send_node(:)

    call setup_sendME(myid,myid_ibegin,myid_npair,myid_npair_orb,loc_bucket(atomf_nlpf_loc))
    call recv_npairME(myid,myid_npair,myid_npair_orb,rem_bucket(atomf_nlpf_loc))
    call setup_recvME(atomf_nlpf_rem,myid,myid_ibegin,rem_bucket(atomf_nlpf_rem),halo(APrange))

    deallocate(isend_array,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating isend_array in set_bucket: ",stat)
    call stop_timer(tmr_std_indexing)
    return
!!***

  contains  ! to the end of this module

!!****f* set_bucket_module/calc_mx_nodes_bucket *
!!
!!  NAME 
!!   calc_mx_nodes_bucket
!!  USAGE
!! 
!!  PURPOSE
!!    calculates the maximum number of nodes for local and remote buckets
!!   In the calculation of matrix elements..., 
!!   Domain Responsible (DR) nodes calculates my domain partial contribution,
!!   stores and sends the values (matrix elements) using their local buckets.
!!   Bundle Responsible (BR) nodes receive the contributions and store them
!!   into my local array using remote buckets.
!!
!!    mx_recv_BRnode: for local bucket ( number of receiving BR nodes in sending) 
!!    mx_send_DRnode: for remote bucket ( number of sending DR nodes in receiving)
!! 
!!  INPUTS
!! 
!!  USES
!! 
!!  AUTHOR
!!   T. Miyazaki
!!  CREATION DATE
!!   19/10/2006
!!  MODIFICATION HISTORY
!!   2016/07/06 17:30 nakata
!!    Renamed comm_in_BG -> comm_in_BtoG and comBG -> comm_naba_blocks_of_atoms
!!
!!  SOURCE
!!
    subroutine calc_mx_nodes_bucket(myid,comm_naba_blocks_of_atoms,mx_recv_BRnode, mx_send_DRnode)

      use datatypes
      use global_module, ONLY: numprocs
      use primary_module, ONLY: bundle
      use naba_blk_module,ONLY: comm_in_BtoG
      use GenComms, ONLY: cq_abort

      implicit none
      integer, intent(in) :: myid
      type(comm_in_BtoG),intent(in) :: comm_naba_blocks_of_atoms
      integer,intent(out) :: mx_send_DRnode, mx_recv_BRnode

      integer :: inp,nnode,iprim,inode,jnode,ind_node,ind_node2,ifind
      integer :: stat,irc,ierr
     !Automatic array
      integer :: list_send_node(numprocs)
      integer :: list_recv_node(numprocs)

      inp=0; nnode=0
      do iprim=1,bundle%mx_iprim
         !do iprim=1,bundle%n_prim
         if(comm_naba_blocks_of_atoms%no_recv_node(iprim) > 0) then
            inp=inp+1
            if(inp == 1) then    ! if this is the first
               nnode=comm_naba_blocks_of_atoms%no_recv_node(iprim)
               do inode=1,comm_naba_blocks_of_atoms%no_recv_node(iprim)
                  list_send_node(inode)= &
                       comm_naba_blocks_of_atoms%list_recv_node(inode,iprim)
               enddo
            else
               do inode=1,comm_naba_blocks_of_atoms%no_recv_node(iprim)
                  ! naba nodes for iprim-th atom
                  ind_node=comm_naba_blocks_of_atoms%list_recv_node(inode,iprim)
                  ifind=0
                  do jnode=1,nnode     ! check existed naba nodes for my bundle
                     ind_node2=list_send_node(jnode)
                     if(ind_node2 == ind_node) then
                        ifind=1
                        exit
                     endif
                  enddo
                  if(ifind == 0) then  ! if this node is new
                     nnode=nnode+1
                     list_send_node(nnode)= comm_naba_blocks_of_atoms%list_recv_node(inode,iprim)
                  endif ! if this node is new
               enddo  ! end do loop over neighbour nodes for iprim-th atom
            endif    ! present primary atom is first or not
         endif     ! Does iprim-th atom have remote neighbour blocks?
      enddo
      mx_send_DRnode=nnode

  !List of receiving BR nodes
  ! In the blip-grid transforms, 
      inp=0; nnode=0
      do iprim=1,bundle%mx_iprim
         !do iprim=1,bundle%n_prim
         if(comm_naba_blocks_of_atoms%no_send_node(iprim) > 0) then
            inp=inp+1
            if(inp == 1) then    ! if this is the first
               nnode=comm_naba_blocks_of_atoms%no_send_node(iprim)
               do inode=1,comm_naba_blocks_of_atoms%no_send_node(iprim)
                  list_recv_node(inode)= &
                       comm_naba_blocks_of_atoms%list_send_node(inode,iprim)
               enddo
            else
               do inode=1,comm_naba_blocks_of_atoms%no_send_node(iprim)
                  ! naba nodes for iprim-th atom
                  ind_node=comm_naba_blocks_of_atoms%list_send_node(inode,iprim)
                  ifind=0
                  do jnode=1,nnode     ! check existed naba nodes for my bundle
                     ind_node2=list_recv_node(jnode)
                     if(ind_node2 == ind_node) then
                        ifind=1
                        exit
                     endif
                  enddo
                  if(ifind == 0) then  ! if this node is new
                     nnode=nnode+1
                     list_recv_node(nnode)= comm_naba_blocks_of_atoms%list_send_node(inode,iprim)
                  endif ! if this node is new
               enddo  ! end do loop over neighbour nodes for iprim-th atom
            endif    ! present primary atom is first or not
         endif     ! Does iprim-th atom have remote neighbour blocks?
      enddo
      mx_recv_BRnode=nnode

      if(iprint_index >= 4.AND.myid==0) write(io_lun,111) myid+1, mx_recv_BRnode, mx_send_DRnode
      111 format(' Calc_mx_send for local and remote buckets : Node = ',i4,&
                 ' mx_recv_BRnode & mx_send_DRnode = ',2i4)
     return
    end subroutine calc_mx_nodes_bucket

!!***

!!****f* set_bucket_module/make_sendinfo_ME *
!!
!!  NAME 
!!   make_sendinfo_ME
!!  USAGE
!! 
!!  PURPOSE
!!    calculates rembucket%no_send_node & 
!!               rembucket%list_send_node(mx_send_node)
!!    to prepare receiving MPI calls
!! 
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   T. Miyazaki
!!  CREATION DATE
!!   30/05/2000
!!  MODIFICATION HISTORY
!!   2006/07/03 08:05 dave
!!    Various changes associated with variable NSF and general tidying
!!   2016/07/06 17:30 nakata
!!    Renamed comm_in_BG -> comm_in_BtoG and comBG -> comm_naba_blocks_of_atoms
!!  SOURCE
!!
    subroutine make_sendinfo_ME(rembucket,comm_naba_blocks_of_atoms)

      use datatypes
      use primary_module, ONLY: bundle
      use bucket_module,  ONLY: remote_bucket  ! Can we write in this way?
      use naba_blk_module,ONLY: comm_in_BtoG
      use GenComms, ONLY: cq_abort

      implicit none
      type(comm_in_BtoG),intent(in) :: comm_naba_blocks_of_atoms
      type(remote_bucket),intent(inout):: rembucket

      integer :: inp,nnode,iprim,inode,jnode,ind_node,ind_node2,ifind
      integer :: stat,irc,ierr

      inp=0
      do iprim=1,bundle%mx_iprim
         !do iprim=1,bundle%n_prim
         if(comm_naba_blocks_of_atoms%no_recv_node(iprim) > 0) then
            inp=inp+1
            if(inp == 1) then    ! if this is the first 
               nnode=comm_naba_blocks_of_atoms%no_recv_node(iprim)
               rembucket%no_send_node=nnode
               do inode=1,comm_naba_blocks_of_atoms%no_recv_node(iprim)
                  rembucket%list_send_node(inode)= &
                       comm_naba_blocks_of_atoms%list_recv_node(inode,iprim)
               enddo
            else
               do inode=1,comm_naba_blocks_of_atoms%no_recv_node(iprim) 
                  ! naba nodes for iprim-th atom
                  ind_node=comm_naba_blocks_of_atoms%list_recv_node(inode,iprim)
                  ifind=0
                  do jnode=1,nnode     ! check existed naba nodes for my bundle
                     ind_node2=rembucket%list_send_node(jnode)
                     if(ind_node2 == ind_node) then
                        ifind=1
                        exit
                     endif
                  enddo
                  if(ifind == 0) then  ! if this node is new
                     rembucket%no_send_node=rembucket%no_send_node+1
                     rembucket%list_send_node(rembucket%no_send_node)= &
                          comm_naba_blocks_of_atoms%list_recv_node(inode,iprim)
                     nnode=rembucket%no_send_node
                  endif ! if this node is new
               enddo  ! end do loop over neighbour nodes for iprim-th atom
            endif    ! present primary atom is first or not
         endif     ! Does iprim-th atom have remote neighbour blocks?
      enddo

      !check remote_buc
      if(rembucket%no_send_node > rembucket%mx_send_node) &
           call cq_abort('ERROR! in make_sendinfo_ME: no_send_nodes = ', &
           & rembucket%no_send_node,rembucket%mx_send_node)
      return
    end subroutine make_sendinfo_ME
!!***

!!****f* set_bucket_module/make_pair_DCSpart *
!!
!!  NAME 
!!   make_pair_DCSpart
!!  USAGE
!! 
!!  PURPOSE
!!    makes the information of pairs of two atoms, whose distances
!!   are within a certain cutoff radius.
!!    For these pairs, domain partial contributions of matrix elements
!!   will be stored and sent to their responsible bundle responsible
!!   nodes. (local_bucket)
!! 
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   T. Miyazaki
!!  CREATION DATE
!!   30/05/2000
!!  MODIFICATION HISTORY
!!   2006/07/03 08:40 dave
!!    Various variable NSF changes and tidying
!!   2016/07/20 16:30 nakata
!!    Renamed naba_atm -> naba_atoms_of_blocks, halo_atm -> halo_atoms_of_blocks
!!  SOURCE
!!
    subroutine make_pair_DCSpart(loc_bucket,max_pair)

      use datatypes
      use numbers,        ONLY: RD_ERR
      use global_module,  ONLY: numprocs, id_glob, species_glob, sf, nlpf, paof
      use species_module, ONLY: nsf_species, nlpf_species, npao_species
      use group_module,   ONLY: parts
      use primary_module, ONLY: domain
      use cover_module,   ONLY: DCS_parts,BCS_parts
      use naba_blk_module,ONLY: naba_atm_of_blk,halo_atm_of_blk
      use set_blipgrid_module, ONLY: naba_atoms_of_blocks, halo_atoms_of_blocks
      use GenComms, ONLY: cq_abort, myid

      implicit none

      type(local_bucket),intent(inout) :: loc_bucket
      integer, intent(out) :: max_pair

      type(naba_atm_of_blk),pointer :: naba_atm1,naba_atm2
      type(halo_atm_of_blk),pointer :: halo_atm1,halo_atm2

      real(double) :: rcutsq,distsq
      real(double) :: xmu,ymu,zmu,xnu,ynu,znu
      integer :: iprim_blk,mpart,inp1,np1,mu,ni1,ind_halo1
      integer :: inp2,npart,np2,nu,ni2,ind_halo2,npair,ind_pair,norb, ind_orb,norb1,norb2
      integer :: ind_cover1,ind_cover2,icover1,icover2,ip1,ip2
      integer :: nnd_old1,no_of_node,nnd_rem,nnd_old
      integer :: ncoverz,ncoveryz,ncoverz_rem,ncoveryz_rem,ind_part
      integer :: nx,ny,nz,offset,i
      integer :: stat,irc,ierr
      integer :: nx1,nx2,ny1,ny2,nz1,nz2, spec1, spec2

      npair=0
      !TM VARNSF: START
      norb = 0
      !TM VARNSF: END
      rcutsq=loc_bucket%rcut*loc_bucket%rcut
      if(rcutsq < RD_ERR) call cq_abort('ERROR in make_pair_DCSpart : rcutsq= ',rcutsq)

      !local bucket should identify left and right functions
      !-- Now the following two lines have already been done 
      !   in set_local_bucket. (14/Jul/2000 T. Miyazaki)
      !OLD  loc_bucket%naba_atm_left => naba_atm1
      !OLD  loc_bucket%naba_atm_right => naba_atm2
      naba_atm1 => naba_atoms_of_blocks(loc_bucket%ind_left)
      naba_atm2 => naba_atoms_of_blocks(loc_bucket%ind_right)
      halo_atm1 => halo_atoms_of_blocks(loc_bucket%ind_left)
      halo_atm2 => halo_atoms_of_blocks(loc_bucket%ind_right)

      !initialize table
      loc_bucket%i_h2d=0

      !making table of pairs of two neighbour atoms for all primary blocks
      ! just check which pairs(halo1,halo2) are neigbour pairs.
      do iprim_blk=1,domain%groups_on_node
         if(naba_atm1%no_of_part(iprim_blk) > 0) then
            inp1=0
            if(iprim_blk > naba_atm1%mx_iprim_blk) &
                 call cq_abort('iprim_blk error1 in make_pair ',iprim_blk,naba_atm1%mx_iprim_blk)
            do mpart=1,naba_atm1%no_of_part(iprim_blk)    !naba parts (NOPG)
               np1=naba_atm1%list_part(mpart,iprim_blk)   !NOPG label in cov. set
               if(mpart > naba_atm1%mx_part) &
                    call cq_abort('mpart error in make_pair ',mpart,naba_atm1%mx_part)
               do mu=1,naba_atm1%no_atom_on_part(mpart,iprim_blk)!naba atms in part
                  inp1=inp1+1
                  if(inp1 > naba_atm1%mx_atom) &
                       call cq_abort('inp1 error in make_pair ',inp1,naba_atm1%mx_atom)
                  ni1=naba_atm1%list_atom(inp1,iprim_blk)
                  ind_halo1=naba_atm1%list_atom_by_halo(inp1,iprim_blk)

                  if(np1 > DCS_parts%ng_cover ) &
                       call cq_abort(' np1 error in make_pair ',np1,DCS_parts%ng_cover )
                  if(ni1 > DCS_parts%n_ing_cover(np1) ) &
                       call cq_abort(' ni1 error in make_pair ',ni1,DCS_parts%n_ing_cover(np1))
                  if(ind_halo1 > loc_bucket%no_halo_atom1) &
                       call cq_abort('ind_halo1 error ',ind_halo1,loc_bucket%no_halo_atom1)
                  if(ind_halo1 < 1 ) &
                       call cq_abort('ind_halo1 is not halo atoms!! ',inp1,ind_halo1)
                  xmu=DCS_parts%xcover(DCS_parts%icover_ibeg(np1)+ni1-1)
                  ymu=DCS_parts%ycover(DCS_parts%icover_ibeg(np1)+ni1-1)
                  zmu=DCS_parts%zcover(DCS_parts%icover_ibeg(np1)+ni1-1)

                  if(naba_atm2%no_of_part(iprim_blk) > 0) then
                     inp2=0
                     do npart=1,naba_atm2%no_of_part(iprim_blk)  ! naba parts for right
                        np2=naba_atm2%list_part(npart,iprim_blk) ! NOPG lab in cov. set
                        if(npart > naba_atm2%mx_part) &
                             call cq_abort('npart error in make_pair ',npart,naba_atm2%mx_part)
                        do nu=1,naba_atm2%no_atom_on_part(npart,iprim_blk) ! naba atms 
                           inp2=inp2+1
                           if(inp2 > naba_atm2%mx_atom) &
                                call cq_abort('inp2 error in make_pair ',inp2,naba_atm2%mx_atom)
                           ni2=naba_atm2%list_atom(inp2,iprim_blk)
                           ind_halo2=naba_atm2%list_atom_by_halo(inp2,iprim_blk)

                           if(np2 > DCS_parts%ng_cover ) &
                                call cq_abort(' np2 error in make_pair ',np2,DCS_parts%ng_cover )
                           if(ni2 > DCS_parts%n_ing_cover(np2) ) &
                                call cq_abort(' ni2 error in make_pair ',ni2,DCS_parts%n_ing_cover(np2) )
                           if(ind_halo2 > loc_bucket%no_halo_atom2) &
                                call cq_abort('ind_halo2 error ',ind_halo2,loc_bucket%no_halo_atom2)
                           xnu=DCS_parts%xcover(DCS_parts%icover_ibeg(np2)+ni2-1)
                           ynu=DCS_parts%ycover(DCS_parts%icover_ibeg(np2)+ni2-1)
                           znu=DCS_parts%zcover(DCS_parts%icover_ibeg(np2)+ni2-1)

                        !(VARIOUS-CUTOFF)
                        ! In the future, when we introduce various cutoffs depending on the species of atoms,
                        ! rcutsq in the following needs to be changed by (rcut(i) + rcut(j)) **2
                        !        commented by TM 18Oct2006
   
                           distsq=(xmu-xnu)**2+(ymu-ynu)**2+(zmu-znu)**2
                           if(distsq < rcutsq ) then
                              if(loc_bucket%i_h2d(ind_halo2,ind_halo1) == 0) then !new pair
                                 npair=npair+1
                                 loc_bucket%i_h2d(ind_halo2,ind_halo1)=npair
                                 !TM VARNSF: START
                                 !!     Here, we do not need to change the definition of i_h2d, as
                                 !!     we will redefine this later.

                                 !!   I suggest to have function_type as a member of naba_atm_of_blk
                                 !func_type1 = naba_atm1%function_type
                                 !func_type2 = naba_atm2%function_type
                                 !norb1: number of orbitals for halo atom 1, and func_type1
                                 !norb2: number of orbitals for halo atom 2, and func_type2
                                 ! DRB, 2006/06/28 15:04: func_type is set by the local bucket; 
                                 ! we can get norb1 and norb2 from DCS_parts and lab_cell
                                 spec1 = species_glob( id_glob( parts%icell_beg(DCS_parts%lab_cell(np1)) +ni1-1 ))
                                 spec2 = species_glob( id_glob( parts%icell_beg(DCS_parts%lab_cell(np2)) +ni2-1 ))
                                 select case(loc_bucket%ind_left)
                                 case(sf)
                                    norb1 = nsf_species(spec1)
                                 case(nlpf)
                                    norb1 = nlpf_species(spec1)
                                 case(paof)
                                    norb1 = npao_species(spec1)
                                 end select
                                 select case(loc_bucket%ind_right)
                                 case(sf)
                                    norb2 = nsf_species(spec2)
                                 case(nlpf)
                                    norb2 = nlpf_species(spec2)
                                 case(paof)
                                    norb2 = npao_species(spec2)
                                 end select
                                 norb = norb+norb1*norb2
                                 !TM VARNSF: END
                              endif
                           endif
                        enddo         !loop over nu atoms
                     enddo          !loop over npart partitions
                  endif           ! if there are naba2 parts for this block
               enddo           !loop over mu atoms
            enddo            !loop over mpart partitions
         endif             ! if there are naba1 parts for this block
      enddo             !loop over primary blocks

      ! check npair 
      loc_bucket%no_pair=npair
      !TM VARNSF: START
      loc_bucket%no_pair_orb = norb
      !TM VARNSF: END
      if(iprint_index>3.AND.myid==0) write(io_lun,*) 'n_pair in make_pair_DCSpart = ',npair
      !if(npair > loc_bucket%mx_pair) call cq_abort('ERROR in make_pair_DCSpart for loc_bucket%mx_pair',&
      !     npair,loc_bucket%mx_pair)
      !making tables

!      ind_pair=0
!      ! Count pairs
!      do ip1=1,halo_atm1%no_of_part   ! Loop over halo partitions for left
!         np1=halo_atm1%list_of_part(ip1)        ! NOPG label
!         if(np1 > DCS_parts%ng_cover) call cq_abort(' ERROR :np1 in make_pair_DCSpart ',np1)
!         do ni1=1,DCS_parts%n_ing_cover(np1)     ! do loop over atoms
!            icover1=DCS_parts%icover_ibeg(np1)+ni1-1
!            ind_halo1=halo_atm1%ihalo(icover1)
!            if(ind_halo1 /= 0) then            ! if ni1 is halo atom, then
!               icover2=0
!               do ip2=1,halo_atm2%no_of_part     ! Loop over halo parts for right
!                  np2=halo_atm2%list_of_part(ip2)     !NOPG label in a cov. set
!                  do ni2=1,DCS_parts%n_ing_cover(np2) ! do loop over atoms
!                     icover2=DCS_parts%icover_ibeg(np2)+ni2-1
!                     ind_halo2=halo_atm2%ihalo(icover2)
!                     if(ind_halo2 /= 0) then
!                        if(loc_bucket%i_h2d(ind_halo2,ind_halo1) /= 0) ind_pair=ind_pair+1
!                     end if
!                  end do
!               end do
!            end if
!         end do
!      end do
      ! End counting pairs
      if(allocated(isend_array).AND.myid==0) write(io_lun,*) "isend_array allocated: ",size(isend_array)
      allocate(isend_array(4*npair),STAT=stat)
      if(stat/=0) call cq_abort("Error allocating isend_array in make_pair_DCSprt: ",npair)
      ind_pair=0
      !TM VARNSF: START
      ind_orb=0
      !TM VARNSF: END
      icover1=0
      nnd_old=-1
      no_of_node=0
      ncoverz=DCS_parts%ncoverz
      ncoveryz=DCS_parts%ncovery*DCS_parts%ncoverz
      if(halo_atm1%no_of_part > 0) then
         do ip1=1,halo_atm1%no_of_part   ! Loop over halo partitions for left
            np1=halo_atm1%list_of_part(ip1)        ! NOPG label
            if(np1 > DCS_parts%ng_cover) call cq_abort(' ERROR :np1 in make_pair_DCSpart ',np1)
            ind_part=DCS_parts%lab_cell(np1)       ! CClabel in a sim. cell
            !OLD if(ind_part > DCS_parts%ng_cover) then
            if(ind_part > parts%mx_gcell) call cq_abort(' ERROR :ind_part in make_pair_DCSpart ',ind_part, DCS_parts%ng_cover)
            ind_cover1=DCS_parts%lab_cover(np1)    ! CClabel in a cov. set
            nx1=(ind_cover1-1)/(DCS_parts%ncovery*DCS_parts%ncoverz)
            ny1=(ind_cover1-1-nx1*DCS_parts%ncovery*DCS_parts%ncoverz)&
                 /DCS_parts%ncoverz
            nz1=ind_cover1-nx1*DCS_parts%ncovery*DCS_parts%ncoverz&
                 -ny1*DCS_parts%ncoverz

            if(ind_cover1 > DCS_parts%ng_cover) call cq_abort(' ERROR :ind_cover1 in make_pair_DCSpart ',ind_cover1)
            nnd_rem=parts%i_cc2node(ind_part)
            if(nnd_rem > numprocs) call cq_abort(' ERROR :nnd_rem in make_pair_DCSpart ',nnd_rem)

            !IF new node then
            if(nnd_rem /= nnd_old) then
               no_of_node=no_of_node+1
               if(no_of_node > loc_bucket%mx_recv_node) &
                    call cq_abort('ERROR in make_pair_DCSpart for loc_bucket rec nodes', no_of_node, loc_bucket%mx_recv_node)

               loc_bucket%ibegin_pair(no_of_node)=ind_pair+1
               !TM VARNSF: START
               loc_bucket%ibegin_pair_orb(no_of_node)=ind_orb+1
               !TM VARNSF: END
               loc_bucket%list_recv_node(no_of_node)=nnd_rem
               ncoverz_rem =BCS_parts%ncover_rem(3+(nnd_rem-1)*3)
               ncoveryz_rem=BCS_parts%ncover_rem(2+(nnd_rem-1)*3)* &
                    BCS_parts%ncover_rem(3+(nnd_rem-1)*3)
            endif
            nnd_old=nnd_rem
            !endif for new node case

            do ni1=1,DCS_parts%n_ing_cover(np1)     ! do loop over atoms
               icover1=DCS_parts%icover_ibeg(np1)+ni1-1
               ind_halo1=halo_atm1%ihalo(icover1)
               if(ind_halo1 /= 0) then            ! if ni1 is halo atom, then
                  icover2=0
                  if(halo_atm2%no_of_part > 0) then
                     do ip2=1,halo_atm2%no_of_part     ! Loop over halo parts for right
                        np2=halo_atm2%list_of_part(ip2)     !NOPG label in a cov. set
                        ind_cover2=DCS_parts%lab_cover(np2) !CC label in a cov. set
                        nx2=(ind_cover2-1)/(DCS_parts%ncovery*DCS_parts%ncoverz)
                        ny2=(ind_cover2-1-nx2*DCS_parts%ncovery*DCS_parts%ncoverz)&
                             /DCS_parts%ncoverz
                        nz2=ind_cover2-nx2*DCS_parts%ncovery*DCS_parts%ncoverz&
                             -ny2*DCS_parts%ncoverz

                        nx=nx2-nx1
                        ny=ny2-ny1
                        nz=nz2-nz1

                        offset=nx*ncoveryz_rem+ny*ncoverz_rem+nz !offset in remote's cov. set

                        do ni2=1,DCS_parts%n_ing_cover(np2) ! do loop over atoms
                           icover2=DCS_parts%icover_ibeg(np2)+ni2-1
                           ind_halo2=halo_atm2%ihalo(icover2)
                           if(ind_halo2 /= 0) then       ! if ni2 is halo atom, then
                              !the following is executed only if both atoms are halo atoms
                              if(loc_bucket%i_h2d(ind_halo2,ind_halo1) /= 0) then
                                 ! the following is executed only if two atoms make a pair
                                 ! whose domain partial contribution will be calculated
                                 ind_pair=ind_pair+1
                                 if(ind_pair>npair) &
                                      call cq_abort(' pair error: ',ind_pair,npair)
                                 loc_bucket%i_h2d(ind_halo2,ind_halo1)=ind_pair
                                 !TM VARNSF: START
                                 !func_type1 = naba_atm1%function_type
                                 !func_type2 = naba_atm2%function_type
                                 !norb1 <= ind_halo1 (,func_type1)
                                 !norb2 <= ind_halo2, func_type2
                                 ! NOTE !!!!
                                 ! the information of halo atoms must have been made
                                 ! by the specific cutoff, AND the type of functions.
                                 spec1 = species_glob( id_glob( parts%icell_beg(DCS_parts%lab_cell(np1)) +ni1-1 ))
                                 spec2 = species_glob( id_glob( parts%icell_beg(DCS_parts%lab_cell(np2)) +ni2-1 ))
                                 select case(loc_bucket%ind_left)
                                 case(sf)
                                    norb1 = nsf_species(spec1)
                                 case(nlpf)
                                    norb1 = nlpf_species(spec1)
                                 case(paof)
                                    norb1 = npao_species(spec1)
                                 end select
                                 select case(loc_bucket%ind_right)
                                 case(sf)
                                    norb2 = nsf_species(spec2)
                                 case(nlpf)
                                    norb2 = nlpf_species(spec2)
                                 case(paof)
                                    norb2 = npao_species(spec2)
                                 end select
                                 loc_bucket%i_h2d(ind_halo2,ind_halo1)=ind_orb+1
                                 ind_orb =ind_orb + norb1*norb2
                                 !TM VARNSF: END
                                 if(4*ind_pair>4*npair) call cq_abort('make_pair ERROR: ',ind_pair,npair)
                                 isend_array(4*ind_pair-3)=ind_part
                                 isend_array(4*ind_pair-2)=ni1
                                 isend_array(4*ind_pair-1)=ni2
                                 isend_array(4*ind_pair  )=offset
                              endif ! if ni2 is halo atom
                           endif  ! if (ni1,ni2) is a pair
                        enddo ! loop over ni2 (atoms in the halo part)
                     enddo   ! loop over ip2 (halo partitions for right)
                  end if ! (halo_atm2%no_of_part > 0) 
               endif  ! if ni1 is halo_atom
            enddo   ! loop over ni1 (atoms in the halo part)
         enddo    ! loop over ip1 (halo partitions for left)
      end if ! (halo_atm1%no_of_part > 0)
      loc_bucket%no_recv_node = no_of_node
      max_pair = loc_bucket%no_pair
      !TM VARNSF: START
      !I would like to add here the following to check the number of pair of orbitals.
      if(ind_orb /= loc_bucket%no_pair_orb) &
           call cq_abort('ERROR: loc_bucket%no_pair_orb mismatch ',ind_orb, loc_bucket%no_pair_orb)
      if(ind_pair /= loc_bucket%no_pair) &
           call cq_abort('ERROR: loc_bucket%no_pair mismatch ',ind_pair, loc_bucket%no_pair)
      !TM VARNSF: END
      return
    end subroutine make_pair_DCSpart
!!***

!!****f* set_bucket_module/setup_sendME *
!!
!!  NAME 
!!   setup_sendME
!!  USAGE
!! 
!!  PURPOSE
!!    sending the information of pairs of two atoms, for which
!!   domain partial contributions of matrix elements will be
!!   stored and sent to their bundle responsible nodes.
!!    There are two MPI_issend calls in this subroutine. The first one
!!   is just to send the size of the pairs. After having obtained this
!!   information, remote (bundle responsible) nodes can receive
!!   the information of pairs from domain responsible nodes.
!! 
!!   In the calculation of surface etc, we should consider the possibility
!!   of having no receiving nodes.
!!
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!! 
!!  CREATION DATE
!! 
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
    subroutine setup_sendME(myid,myid_ibegin,myid_npair,myid_npair_orb,loc_bucket)

      use datatypes
      use mpi
      use bucket_module, ONLY:local_bucket
      use GenComms, ONLY: cq_abort
      use global_module, ONLY: iprint_index

      implicit none

      integer,intent(in)  :: myid
      !integer,intent(out) :: myid_ibegin,myid_npair
      !TM VARNSF: START
      integer,intent(out) :: myid_ibegin,myid_npair,myid_npair_orb
      !TM VARNSF: END
      type(local_bucket),intent(in) :: loc_bucket

      integer :: mynode,nnd,nnd_rem,ibegin,send_size
      integer :: irc,ierr,tag
      integer :: ii, stat

      mynode=myid+1
      myid_ibegin=0
      myid_npair=0

      ! if(loc_bucket%no_recv_node < 1) then
      !    call cq_abort('ERROR in setup_sendME : no_recv_node = ', &
      !                   mynode, loc_bucket%no_recv_node)
      ! endif
      if(loc_bucket%no_recv_node > 0) then
         do nnd=1,loc_bucket%no_recv_node   ! Loop over receiving nodes
            nnd_rem=loc_bucket%list_recv_node(nnd)
            if(nnd < loc_bucket%no_recv_node) then
               buff_npair(nnd)=&
                    loc_bucket%ibegin_pair(nnd+1)-loc_bucket%ibegin_pair(nnd)
               !TM VARNSF : START
               buff_npair_orb(nnd)=&
                    loc_bucket%ibegin_pair_orb(nnd+1)-loc_bucket%ibegin_pair_orb(nnd)
               !TM VARNSF : END
            else 
               buff_npair(nnd)=loc_bucket%no_pair-loc_bucket%ibegin_pair(nnd)+1
               !TM VARNSF : START
               buff_npair_orb(nnd)=loc_bucket%no_pair_orb-loc_bucket%ibegin_pair_orb(nnd)+1
               !TM VARNSF : END
            endif

            ! keep the pointer (mynode case) or sending arrays (remote node case)
            if(nnd_rem == mynode) then
               myid_ibegin=4*(loc_bucket%ibegin_pair(nnd)-1)+1
               myid_npair=buff_npair(nnd)
               !TM VARNSF : START
               myid_npair_orb=buff_npair_orb(nnd)
                if(iprint_index>4.AND.myid==0) write(io_lun,112) mynode ,myid_ibegin,myid_npair, myid_npair_orb
                112 format(' Node ',i3,' myid_ibegin, myid_npair, myid_npair_orb = ',3i6)
               !if(iprint_index>4) write(io_lun,*) ' Node ',mynode, &
               !     ' myid_ibegin, myid_npair and myid_npair_orb= ',myid_ibegin,myid_npair,myid_npair_orb
               !TM VARNSF : END
               !if(iprint_index>4) write(io_lun,*) ' Node ',mynode, ' isend for my node =',(isend_array(myid_ibegin+ii),ii=0,3)
            else
               !!  first, sending buff_npair to nnd_rem
               !tag=1
               !call MPI_issend(buff_npair(nnd),1,MPI_INTEGER,nnd_rem-1,&
               !     tag,MPI_COMM_WORLD,nsend_req(2*nnd-1),ierr)
               !!  sending loc_bucket%send_array(ibegin:ibegin+send_size-1) to nnd_rem
               !tag=2
               !ibegin=4*(loc_bucket%ibegin_pair(nnd)-1)+1
               !send_size=4*buff_npair(nnd)
               !call MPI_issend(isend_array(ibegin),send_size,MPI_INTEGER,&
               !     nnd_rem-1,tag,MPI_COMM_WORLD,nsend_req(2*nnd),ierr)
               !TM VARNSF : START
               !first, sending buff_npair to nnd_rem
               tag=1
               call MPI_issend(buff_npair(nnd),1,MPI_INTEGER,nnd_rem-1,&
                    tag,MPI_COMM_WORLD,nsend_req(3*nnd-1),ierr)
               !second, sending bunn_npair_orb to nnd_rem (added for variable NSF)
               tag=2
               call MPI_issend(buff_npair_orb(nnd),1,MPI_INTEGER,nnd_rem-1,&
                    tag,MPI_COMM_WORLD,nsend_req(3*nnd-2),ierr)
               !sending loc_bucket%send_array(ibegin:ibegin+send_size-1) to nnd_rem
               tag=3
               ibegin=4*(loc_bucket%ibegin_pair(nnd)-1)+1
               send_size=4*buff_npair(nnd)
               if(send_size>0) then
                  if(ibegin>4*loc_bucket%no_pair.AND.myid==0) write(io_lun,*) 'ERROR: ',ibegin,send_size,4*loc_bucket%no_pair
                  if(ibegin + send_size-1>4*loc_bucket%no_pair.AND.myid==0) &
                       write(io_lun,*) 'ERROR: ',ibegin,send_size,4*loc_bucket%no_pair
                  call MPI_issend(isend_array(ibegin),send_size,MPI_INTEGER,&
                       nnd_rem-1,tag,MPI_COMM_WORLD,nsend_req(3*nnd),ierr)
                  !TM VARNSF : END
               end if
            endif
            !if(iprint_index>4) write(io_lun,*) '$ME$ Node ',myid+1, ' RecvNode ',nnd_rem,' npair= ',buff_npair(nnd)

         enddo ! End do loop over receiving nodes
         !deallocate(buff_npair, buff_npair_orb,STAT=stat)
         !if(stat/=0) call cq_abort('Error deallocating buff_npair in setup_sendME: ', loc_bucket%no_recv_node)
      endif  !(loc_bucket%no_recv_node > 0) 

      return
    end subroutine setup_sendME
!!***

!!****f* set_bucket_module/recv_npairME *
!!
!!  NAME 
!!   recv_npairME
!!  USAGE
!! 
!!  PURPOSE
!!   receives the number of pairs 
!! 
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!! 
!!  CREATION DATE
!! 
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
    subroutine recv_npairME(myid,myid_npair,myid_npair_orb,rem_bucket)

      use mpi, ONLY: MPI_INTEGER, MPI_COMM_WORLD
      use bucket_module, ONLY:remote_bucket, reset_remote_bucket
      use GenComms, ONLY: cq_abort

      implicit none
      !integer,intent(in)::myid,myid_npair
      !TM VARNSF : START
      integer,intent(in)::myid,myid_npair, myid_npair_orb
      integer :: tag1 = 1, tag2=2
      !TM VARNSF : END
      type(remote_bucket),intent(inout) :: rem_bucket

      integer ::mynode,nnd,inode,irc,ierr,tag
      logical :: change

      mynode=myid+1
      tag=1

      if(rem_bucket%no_send_node < 1) call cq_abort('ERROR in recv_npairME : no_send_node=',rem_bucket%no_send_node)
      change = .false.
      do nnd=1,rem_bucket%no_send_node
         inode=rem_bucket%list_send_node(nnd)

         if(inode /= mynode) then
            !call MPI_recv(rem_bucket%no_of_pair(nnd),1,MPI_INTEGER,inode-1,&
            !     tag,MPI_COMM_WORLD,nrecv_stat(:,nnd),ierr)
            !TM VARNSF : START
            !receive no_of_pair
            call MPI_recv(rem_bucket%no_of_pair(nnd),1,MPI_INTEGER,inode-1,&
                 !tag1,MPI_COMM_WORLD,nrecv_stat(:,nnd),ierr)
                 tag1,MPI_COMM_WORLD,nrecv_stat(:,3*nnd-2),ierr)
            !receive no_of_pair_orb
            call MPI_recv(rem_bucket%no_of_pair_orbs(nnd),1,MPI_INTEGER,inode-1,&
                 !tag2,MPI_COMM_WORLD,nrecv_stat(:,nnd),ierr)
                 tag2,MPI_COMM_WORLD,nrecv_stat(:,3*nnd-1),ierr)
            !TM VARNSF : END
         else
            !rem_bucket%no_of_pair(nnd)=myid_npair
            !TM VARNSF : START
            rem_bucket%no_of_pair(nnd)=myid_npair
            rem_bucket%no_of_pair_orbs(nnd)=myid_npair_orb
            !TM VARNSF : END
         endif
         if(rem_bucket%no_of_pair(nnd)>rem_bucket%mx_pair_comm) then
            rem_bucket%mx_pair_comm = rem_bucket%no_of_pair(nnd)
            change = .true.
         end if
      enddo
      if(change) then
         call reset_remote_bucket(rem_bucket)
      end if
      return
    end subroutine recv_npairME
!!***

!!****f* set_bucket_module/setup_recvME *
!!
!!  NAME 
!!   setup_recvME
!!  USAGE
!! 
!!  PURPOSE
!!    receives the information of pairs and 
!!   constructs the remote_bucket%bucket(ipair,inode) which shows
!!   where to put each sent domain partial contribution of 
!!   matrix element into the local array of the matrix.
!!
!!    For the Hamiltonian matrix, local array of the matrix is 
!!   made to store the whole hamiltonian (including non-local parts),
!!   while domain partial contribution of matrix elements are 
!!   calculated only for those of local parts.
!!   Thus, different cutoff radii are used between local and 
!!   remote buckets.
!! 
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   T. Miyazaki
!!  CREATION DATE
!!   30/05/2000
!!  MODIFICATION HISTORY
!!   2006/12/20 08:41 dave
!!    Small tweak for cq_abort
!!  SOURCE
!!
    subroutine setup_recvME(itype,myid,myid_ibegin,rem_bucket1,mat1_halo,rem_bucket2,mat2_halo)

      use mpi
      use datatypes 
      use numbers, ONLY: RD_ERR
      use group_module,ONLY  :parts
      use primary_module,ONLY:bundle
      use cover_module,ONLY  :BCS_parts
      use matrix_module, ONLY:matrix_halo
      use bucket_module
      use GenComms, ONLY: cq_abort
      use global_module, ONLY: iprint_index

      implicit none

      integer,intent(in) :: itype,myid,myid_ibegin
      type(remote_bucket),intent(inout) :: rem_bucket1
      type(remote_bucket),intent(inout),OPTIONAL :: rem_bucket2
      type(matrix_halo),intent(in) :: mat1_halo
      type(matrix_halo),intent(in),OPTIONAL :: mat2_halo

      integer,pointer :: recv_ptr(:)

      integer :: nnd,nnd_send,id_myid,npair,nsize,ipair
      integer :: ind_part,mu,nu,offset,ipart,nprim,icover,jcover,ni2
      integer :: irc,ierr,tag, ibeg_orb1, ibeg_orb2, norb1, norb2
      logical :: ltype
      integer :: ni2_in_halo,jcover_NOPG,ni2_new,nx,ny,nz, stat
      real(double) :: dx,dy,dz
      real(double) :: rcut
      !For MPI_Wait and VarNSF 08Sep2006
      integer :: nwait_stat1(MPI_STATUS_SIZE),nwait_stat2(MPI_STATUS_SIZE)
      integer :: nwait_stat3(MPI_STATUS_SIZE)
      integer :: nnd_rem, itmp

      rcut=rem_bucket1%locbucket%rcut
      id_myid=myid+1

      ltype=(PRESENT(rem_bucket2).AND.PRESENT(mat2_halo))

      !itype is just for checking ltype
      if(itype == 1 .AND. .NOT.ltype) call cq_abort("Type mismatch in setup_recvME: ",itype)
      do nnd=1,rem_bucket1%no_send_node         !Loop over remote nodes
         nnd_send=rem_bucket1%list_send_node(nnd)  !label of the remote node
         npair=rem_bucket1%no_of_pair(nnd)         !no. of pairs from the node
         if(npair>size(rem_bucket1%bucket,1)) then
            itmp = size(rem_bucket1%bucket,1)
            call cq_abort('Overflow in setup_recvME: ',npair,itmp)
         end if
         nsize=npair*4                             !size
         if(nsize>0) then
            allocate(irecv_array(nsize),STAT=stat)
            if(stat/=0) call cq_abort("Error allocating irecv_array in setup_recvME: ",nsize)
            irecv_array = 0
            if(nnd_send == id_myid) then     !if nnd_send is me, no need of MPI
               recv_ptr => isend_array(myid_ibegin:myid_ibegin+nsize-1)
            else                             !receives information from rem. nodes
               tag = 2
               !TM VARNSF : START
               !as we have add MPI_issend for npair_orb, the index of tag for this MPI-call has become 3
               tag=3
               !TM VARNSF : END
               call MPI_recv(irecv_array,nsize,MPI_INTEGER,nnd_send-1, &
                                !tag,MPI_COMM_WORLD,nrecv_stat(:,nnd),ierr)
                    tag,MPI_COMM_WORLD,nrecv_stat(:,3*nnd),ierr)
               if(ierr /= 0) call cq_abort('ERROR in MPI_recv in put_matrix_elements',ierr)
               recv_ptr => irecv_array(1:nsize)
            endif

            if(iprint_index>3.AND.myid==0) write(io_lun,*) '$ME$ Node ',myid+1,' SendNode ',nnd_send,' npair= ',npair
            !TM VARNSF : START
            ibeg_orb1 = 0
            ibeg_orb2 = 0
            !TM VARNSF : START
            do ipair=1,npair
               ind_part=recv_ptr(4*ipair-3)  ! CC in a sim. cell
               mu=recv_ptr(4*ipair-2)  ! partition seq. no.
               nu=recv_ptr(4*ipair-1)  ! partition seq. no.
               offset=recv_ptr(4*ipair  )  ! offset of two partitions
               !  which include two atoms

               ipart=parts%i_cc2seq(ind_part)  !prim seq. label of this partition
               !ORI if(parts%i_cc2node(ind_part) /= myid+1) call cq_abort(' ERROR : ind_part in setup_recvME ', &
               !ORI      parts%i_cc2node(ind_part),ipair)
               if(parts%i_cc2node(ind_part) /= myid+1) then
                  write(io_lun,*) ' ind_part, mu, nu, offset = ',ind_part, mu, nu,offset
                  write(io_lun,*) ' myid, nnd, nnd_send, npair - ',myid, nnd, nnd_send, npair
                  call cq_abort(' ERROR : ind_part in setup_recvME ', parts%i_cc2node(ind_part),ipair)
               endif
               if(mu > bundle%nm_nodgroup(ipart).or.mu<1) call cq_abort(' ERROR : mu in setup_recvME ',mu,ipair)
               if(ipart > bundle%groups_on_node.or.ipart < 1) call cq_abort(' ERROR : ipart in setup_recvME ',ipart,ipair)
               nprim=bundle%nm_nodbeg(ipart)+mu-1   !prim seq. no. of the left atom
               if(nprim > bundle%n_prim.or.nprim < 1) call cq_abort(' ERROR : nprim in setup_recvME ',nprim,ipair)
               !     icover=bundle%iprim_group(nprim)  !CC in cov. set of the left part
               ! DRB, 11/05/01 to reflect updates of basic types
               icover=BCS_parts%iprim_group(nprim)  !CC in cov. set of the left part
               if(icover > BCS_parts%ng_cover.or.icover < 1) call cq_abort(' ERROR : icover in setup_recvME ',icover,ipair)
               if(ind_part /= BCS_parts%lab_cell(BCS_parts%inv_lab_cover(icover))) &
                    call cq_abort(' ERROR : icover & ind_part in setup_recvME ',&
                    BCS_parts%lab_cell(BCS_parts%inv_lab_cover(icover)), ind_part)
               jcover=icover+offset            !CC in cov. set of the right part
               if(jcover > BCS_parts%ng_cover.or.jcover < 1) call cq_abort(' ERROR : jcover in setup_recvME ',jcover,ipair)
               ni2=mat1_halo%i_hbeg(jcover)+nu-1          !
               jcover_NOPG=BCS_parts%inv_lab_cover(jcover)
               ni2_new=BCS_parts%icover_ibeg(jcover_NOPG)+nu-1

               dx=bundle%xprim(nprim)-BCS_parts%xcover(ni2_new)
               dy=bundle%yprim(nprim)-BCS_parts%ycover(ni2_new)
               dz=bundle%zprim(nprim)-BCS_parts%zcover(ni2_new)
               nx=offset/(BCS_parts%ncovery*BCS_parts%ncoverz)
               ny=(offset-nx*(BCS_parts%ncovery*BCS_parts%ncoverz))/BCS_parts%ncoverz
               nz=offset-nx*(BCS_parts%ncovery*BCS_parts%ncoverz)-ny*BCS_parts%ncoverz
               if(dx*dx+dy*dy+dz*dz > rcut*rcut.AND.myid==0) then
                  if(iprint_index>3) write(io_lun,101) myid+1,ipair,dx,dy,dz,bundle%xprim(nprim)&
                       ,bundle%yprim(nprim),bundle%zprim(nprim)&
                       ,BCS_parts%xcover(ni2_new),BCS_parts%ycover(ni2_new)&
                       ,BCS_parts%zcover(ni2_new)&
                       ,offset,nx,ny,nz
101               format('Node,ipair= ',2i3,' diff =',3f8.3,2x,' prim =',3f8.3,2x,&
                       'cover=',3f8.3,2x,'offset=',i4,' nx,y,z=',3i3)
               endif

               if(ni2 /= ni2_new) call cq_abort(' ERROR : ni2 and ni2_new mismatch !! ' ,ni2,ni2_new)
               if(nu > parts%nm_group(BCS_parts%lab_cell(jcover_NOPG)).or. nu < 1) &
                    call cq_abort(' ERROR : nu in setup_recvME ',nu, parts%nm_group(BCS_parts%lab_cell(jcover_NOPG)))
               if(ni2 > BCS_parts%mx_mcover.or. ni2 < 1) call cq_abort(' ERROR : ni2 in setup_recvME ',ni2,ipair)
               ni2_in_halo=mat1_halo%i_halo(ni2)          ! halo-atom_ seq. no.
               if(ni2_in_halo == 0) call cq_abort('ERROR : non-halo atom1 !  in setup_recvME ' ,ni2,nu)
               if(ni2_in_halo > mat1_halo%ni_in_halo) &
                    call cq_abort(' ERROR : halo_seq number of the atom in setup_recvME ', &
                    ni2_in_halo,mat1_halo%ni_in_halo)
               !rem_bucket1%bucket(ipair,nnd) = mat1_halo%i_h2d( (nprim-1)*mat1_halo%ni_in_halo + ni2_in_halo )

               !TM VARNSF : START
               rem_bucket1%bucket(ipair,nnd)%iprim = nprim
               rem_bucket1%bucket(ipair,nnd)%jhalo = ni2_in_halo
               rem_bucket1%bucket(ipair,nnd)%ibegin_pair = ibeg_orb1+1
               norb1 = mat1_halo%ndimi(nprim)
               norb2 = mat1_halo%ndimj(ni2_in_halo)
               ibeg_orb1=ibeg_orb1+norb1*norb2
               !TM VARNSF : END

               if(ltype) then
                  if(ni2 /= mat2_halo%i_hbeg(jcover)+nu-1) call cq_abort(' ERROR in mat1_halo%i_hbeg and mat2_halo%i_hbeg'&
                       ,ni2,mat2_halo%i_hbeg(jcover)+nu-1)
                  ni2_in_halo=mat2_halo%i_halo(ni2)          ! halo-atom_ seq. no.
                  if(ni2_in_halo == 0) call cq_abort(' ERROR : the present atom is not halo atom2?!  in setup_recvME '&
                       ,ni2_in_halo,ipair)
                  if(ni2_in_halo > mat2_halo%ni_in_halo) call cq_abort(' ERROR : halo_seq number of the atom2 in setup_recvME '&
                       ,ni2_in_halo,mat2_halo%ni_in_halo)
                  !rem_bucket2%bucket(ipair,nnd) = mat2_halo%i_h2d(&
                  !     (nprim-1)*mat2_halo%ni_in_halo + ni2_in_halo )
                  !TM VARNSF : START
                  rem_bucket2%bucket(ipair,nnd)%iprim = nprim
                  rem_bucket2%bucket(ipair,nnd)%jhalo = ni2_in_halo
                  rem_bucket2%bucket(ipair,nnd)%ibegin_pair = ibeg_orb2+1
                  norb1 = mat2_halo%ndimi(nprim)
                  norb2 = mat2_halo%ndimj(ni2_in_halo)
                  ibeg_orb2=ibeg_orb2+norb1*norb2
                  !TM VARNSF : END

               endif
            enddo
            deallocate(irecv_array,STAT=stat)
            if(stat/=0) call cq_abort("Error deallocating irecv_array in setup_recvME: ",nsize)
         end if
      enddo
      ! MPI_Wait for MPI_issend in setup_sendME  ----------
      !       30Nov2005   Tsuyoshi Miyazaki
      if(rem_bucket1%locbucket%no_recv_node > 0) then
         do nnd=1, rem_bucket1%locbucket%no_recv_node
            nnd_rem = rem_bucket1%locbucket%list_recv_node(nnd)
            if(nnd < rem_bucket1%locbucket%no_recv_node) then
               npair= rem_bucket1%locbucket%ibegin_pair(nnd+1)-rem_bucket1%locbucket%ibegin_pair(nnd)
            else
               npair= rem_bucket1%locbucket%no_pair-rem_bucket1%locbucket%ibegin_pair(nnd)+1
            endif
            if(nnd_rem /= id_myid) then
               call MPI_Wait(nsend_req(3*nnd-1), nwait_stat1, ierr)
               if(ierr /= 0) call cq_abort(' ERROR in MPI_Wait1 in setup_recvME ',ierr)
               call MPI_Wait(nsend_req(3*nnd-2), nwait_stat2, ierr)
               if(ierr /= 0) call cq_abort(' ERROR in MPI_Wait2 in setup_recvME ',ierr)
               if(npair > 0) then
                  call MPI_Wait(nsend_req(3*nnd), nwait_stat3, ierr)
                  if(ierr /= 0) call cq_abort(' ERROR in MPI_Wait3 in setup_recvME ',ierr)
               endif
            endif  ! (nnd_rem /= id_myid) then
         enddo ! nnd=1, loc_bucket%no_recv_node
      endif ! (loc_bucket%no_recv_node > 0) then
      ! MPI_Wait for MPI_issend in setup_sendME  ----------

      return
    end subroutine setup_recvME
!!***

  end subroutine set_bucket

!!****f* set_bucket_module/free_bucket *
!!
!!  NAME 
!!   free_bucket
!!  USAGE
!! 
!!  PURPOSE
!!   Frees memory
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   08:08, 08/01/2003 dave
!!  MODIFICATION HISTORY
!!   2016/12/29 19:00 nakata
!!    Removed rem_bucket(4) and loc_bucket(3), which are no longer needed
!!
!!  SOURCE
!!
  subroutine free_bucket

    use bucket_module, ONLY: free_remote_bucket, free_local_bucket

    call free_remote_bucket(rem_bucket(3))
    call free_remote_bucket(rem_bucket(2))
    call free_remote_bucket(rem_bucket(1))
    call free_local_bucket(loc_bucket(2))
    call free_local_bucket(loc_bucket(1))
    return
  end subroutine free_bucket
!!***
end module set_bucket_module
