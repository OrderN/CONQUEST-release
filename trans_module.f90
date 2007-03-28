! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id: trans_module.f90,v 1.3.2.1 2006/03/31 13:08:11 drb Exp $
! ------------------------------------------------------------------------------
! Module trans_module
! ------------------------------------------------------------------------------
! Code area 2: matrices
! ------------------------------------------------------------------------------

!!****h* Conquest/trans_module *
!!  NAME
!!   trans_module
!!  PURPOSE
!!   Routines associated with taking the global transpose of matrices
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   24/01/00
!!  MODIFICATION HISTORY
!!   03/03/00 by D.R.Bowler
!!    Finishing up transposes
!!   12/04/00 by DRB 
!!    Put itran_addr into trans_remote type
!!   20/04/00 by DRB 
!!    Changed send_pair_info to not use pairs type
!!   21/06/2001 dave
!!    Added ROBODoc header and RCS Id and Log tags and cq_abort
!!   24/06/2002 dave
!!    Changed index_transpose to give exact array size and added RCS static tag
!!***
module trans_module

  use GenComms, ONLY: cq_abort, my_barrier

  implicit none

  ! -------------------------------------------------------
  ! RCS ident string for object file id
  ! -------------------------------------------------------
  character(len=80), private :: RCSid = "$Id: trans_module.f90,v 1.3.2.1 2006/03/31 13:08:11 drb Exp $"

contains

!!****f* trans_module/a_and_b *
!!
!!  NAME 
!!   a_and_b
!!  USAGE
!! 
!!  PURPOSE
!!   makes tables to identify a and b: a is the halo atom on the 
!!   sending node (hence the primary atom on the receiving node), 
!!   and b is the primary atom on the sending node (hence the halo 
!!   atom on the receiving node).
!!   For each primary - halo pair, the tables are (within adata):
!!     ipart_a():   sim-cell CC label of a partition
!!     iseq_a():    a seq. no. in his partition
!!     ipart_b(): CC label of b partition relative to a
!!                    partition in GCS of receiving node
!!     iseq_b():  b seq. no. in her partition
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   24/01/00
!!  MODIFICATION HISTORY
!!   21/06/2001 dave
!!    Added ROBODoc header
!!  SOURCE
!!
  subroutine a_and_b(mynode,trans,halo,data,pairs,prim,gcs)

    use datatypes
    use matrix_module
    use basic_types

    implicit none

    ! Passed variables
    type(matrix_trans) :: trans 
    type(matrix_halo)  :: halo
    type(trans_remote) :: data
    type(pair_data)    :: pairs(:)
    type(primary_set)  :: prim
    type(cover_set)    :: gcs

    integer(integ) :: mynode

    ! Local variables
    integer(integ) :: na,ia,npr,nr,nnd_rem,ncoveryz,ncoveryz_rem
    integer(integ) :: np, ind_part,ind_covera,nxa,nya,nza,nxb,nyb,nzb
    integer(integ) :: n,ind_seq,i,ip,ind_coverb,nxd,nyd,nzd,nnpr

    ! Make a useful constant
    ncoveryz=gcs%ncovery*gcs%ncoverz
    na=0
    ia=0
    npr=0
    ! Make the _pair,_seq arrays
    data%i_nd_pair_addr(1)=1
    do nr=1,data%n_rem_node  ! Loop over remote nodes
       pairs(nr)%length = 0
       nnpr=0
       nnd_rem=data%list_rem_node(nr)
       ncoveryz_rem=gcs%ncover_rem(2+(nnd_rem-1)*3)* &
            gcs%ncover_rem(3+(nnd_rem-1)*3)
       do np=1,data%nhp_for_node(nr) ! Loop over halo partitions for this node
          na=na+1
          ind_part=halo%lab_hcell(na)
          ind_covera=halo%lab_hcover(na)-1
          nxa=ind_covera/ncoveryz
          nya=(ind_covera-nxa*ncoveryz)/gcs%ncoverz
          nza=ind_covera-nxa*ncoveryz-nya*gcs%ncoverz
          do n=1,halo%nh_part(na)  ! Loop over atoms in partitions
             ia=ia+1 ! Goes over all atoms in halo 
             ind_seq=halo%j_seq(ia)
             do i=1,trans%n_hnab(ia) ! Loop over neighbours of each atom
                npr=npr+1
                nnpr=nnpr+1
                pairs(nr)%ipart_a(nnpr)=ind_part  ! Local halo atom, a
                pairs(nr)%iseq_a(nnpr)=ind_seq
                ip=trans%i_prim(npr)
                pairs(nr)%submat(nnpr)=halo%ndimj(ia)*halo%ndimi(ip)
                pairs(nr)%length = pairs(nr)%length + pairs(nr)%submat(nnpr)
                ! iprim_part gives the GCS partn (CC-label) of primary atom b
                ind_coverb=gcs%iprim_group(ip)-1  ! Local primary atom, b
                nxb=ind_coverb/ncoveryz
                nyb=(ind_coverb-nxb*ncoveryz)/gcs%ncoverz
                nzb=ind_coverb-nxb*ncoveryz-nyb*gcs%ncoverz
                ! nxd is the offset of the partition B from partition A
                nxd=nxb-nxa
                nyd=nyb-nya
                nzd=nzb-nza
                pairs(nr)%ipart_b(nnpr)=nxd*ncoveryz_rem+ &
                     nyd*gcs%ncover_rem(3+(nnd_rem-1)*3)+nzd
                ! iprim_seq gives the seq number in its partn of primary atom b
                pairs(nr)%iseq_b(nnpr)=prim%iprim_seq(ip)
             enddo
          enddo
       enddo
      if(nr.lt.data%n_rem_node) then
         data%i_nd_pair_addr(nr+1)=data%i_nd_pair_addr(nr)+pairs(nr)%length
      endif
    enddo
    return
  end subroutine a_and_b
!!***

!!****f* trans_module/send_pair_info *
!!
!!  NAME 
!!   send_pair_info
!!  USAGE
!! 
!!  PURPOSE
!!   Send arrays needed by sbrt index_transpose. Data to be sent 
!!   are the tables ipart_a(), iseq_a(), ipart_b(), and iseq_b(), 
!!   contained within the data structure. For the meanings of these, 
!!   see sbrts a_and_b and index_transpose.
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!! 
!!  MODIFICATION HISTORY
!!   21/06/2001 dave
!!    Added ROBODoc header and cq_abort
!!  SOURCE
!!
  subroutine send_pair_info(mynode,pairsind,trans)

    use mpi
    use datatypes
    use matrix_module
    use matrix_comms_module, ONLY : mx_msg_per_part

    implicit none

    ! Passed variables
    type(trans_remote) :: trans
    integer(integ) :: pairsind(:)
    integer(integ) :: mynode

    ! Local variables
    integer(integ) :: nr,nnd_rem,tag,nreq,ierr,nn,irc

    ! Send pairs data to halo nodes
    do nr=1,trans%n_rem_node
       nnd_rem=trans%list_rem_node(nr)
       if(nnd_rem/=mynode) then
          tag=(nnd_rem-1)*mx_msg_per_part
          call MPI_issend(pairsind(4*(trans%i_pair_addr(nr)-1)+1), &
               4*trans%n_pair(nr), &
               MPI_INTEGER,nnd_rem-1,tag+1,MPI_COMM_WORLD,nreq,ierr)
          if(ierr/=MPI_Success) call cq_abort('send_pair_info: error sending')
          call MPI_issend(trans%i_nd_pair_addr(nr),1,MPI_INTEGER, &
               nnd_rem-1,tag+2,MPI_COMM_WORLD,nreq,ierr)
          if(ierr/=MPI_Success) call cq_abort('send_pair_info: error sending')
       endif
    enddo
    call my_barrier ! Make sure everyone's here
    return
  end subroutine send_pair_info
!!***

!!****f* trans_module/index_transpose *
!!
!!  NAME 
!!   index_transpose
!!  USAGE
!! 
!!  PURPOSE
!!   Fetches the remote-address table i_pair_rem() and creates 
!!   the address table itran_addr() used in sbrt mat_tran.
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!! 
!!  MODIFICATION HISTORY
!!   21/06/2001 dave
!!    Added ROBODoc header and cq_abort
!!   24/06/2002 dave
!!    Changed copy of ind to give array dimensions explicitly
!!   2006/09/07 11:51 dave
!!    Added halo_rem to use the halo of the transposed matrix
!!  SOURCE
!!
  subroutine index_transpose(mx_iprim,mynode,data,halo, &
       halo_rem,pairsind,pairs,parts,prim,gcs)

    use mpi
    use matrix_module
    use basic_types
    use matrix_comms_module, ONLY : mx_msg_per_part

    implicit none

    ! Passed variables
    integer(integ) :: mx_iprim,mynode
    integer(integ) :: pairsind(:)
    type(pair_data)    :: pairs(:)
    type(trans_remote) :: data
    type(matrix_halo)  :: halo, halo_rem
    type(group_set) :: parts
    type(primary_set) :: prim
    type(cover_set) :: gcs


    ! Local variables
    integer(integ) :: nnd_rem,tag,nn,ic,nr,iiip,ip,ierr
    integer(integ) :: noffx,noffy,noffz,jjjp,ih,stat
    integer :: nrstat(MPI_STATUS_SIZE)
    ! Temporary local storage
    integer(integ), target :: ind(4*data%mx_pair)
    integer(integ), pointer :: ipart_a_rem(:),iseq_a_rem(:)
    integer(integ), pointer :: ipart_b_rem(:),iseq_b_rem(:)

    ind = 0
    do nr=1,data%n_rem_node  ! Loop over relevant remote nodes
       ipart_a_rem => ind(1:data%n_pair(nr))
       iseq_a_rem  => ind(data%n_pair(nr)+1:2*data%n_pair(nr))
       ipart_b_rem => ind(2*data%n_pair(nr)+1:3*data%n_pair(nr))
       iseq_b_rem  => ind(3*data%n_pair(nr)+1:4*data%n_pair(nr))
       nnd_rem=data%list_rem_node(nr)
       ! Receive info about atom pairs 
       ipart_a_rem = 0
       ipart_b_rem = 0
       iseq_a_rem = 0
       iseq_b_rem = 0
       tag = (mynode-1)*mx_msg_per_part
       if(nnd_rem/=mynode) then
          call MPI_recv(ind,4*data%n_pair(nr),MPI_INTEGER,&
               nnd_rem-1,tag+1,MPI_COMM_WORLD,nrstat,ierr)
          if(ierr/=MPI_Success) call cq_abort('index_trans: error in MPI_recv')
          call MPI_recv(data%rem_pair_addr(nr),1,MPI_INTEGER, &
               nnd_rem-1,tag+2,MPI_COMM_WORLD,nrstat,ierr)
          if(ierr/=MPI_Success) call cq_abort('index_trans: error in MPI_recv')
       else
          ind(1:4*data%n_pair(nr)) = pairsind(4*(data%i_pair_addr(nr)-1)+1: 4*(data%i_pair_addr(nr)-1)+4*data%n_pair(nr))
       endif
       ! Make address itran_addr() for start of data for each pair 
       ic=data%i_pair_addr(nr)
       do nn=1,data%n_pair(nr)
          iiip=parts%i_cc2seq(ipart_a_rem(nn))
          ip=prim%nm_nodbeg(iiip)+iseq_a_rem(nn)-1
          noffx=1+gcs%nspanlx+prim%idisp_primx(iiip)
          noffy=1+gcs%nspanly+prim%idisp_primy(iiip)
          noffz=1+gcs%nspanlz+prim%idisp_primz(iiip)
          ! NB jjjp = jpart within m_kern_max
          jjjp=ipart_b_rem(nn)+ &
               (noffx-1)*gcs%ncovery*gcs%ncoverz+(noffy-1)*gcs%ncoverz+noffz
          if(jjjp>gcs%mx_gcover.OR.jjjp<1) then
             call cq_abort('index_trans: jjjp wrong ',jjjp,gcs%mx_gcover)
          endif
          if(halo_rem%i_hbeg(jjjp)+iseq_b_rem(nn)-1>gcs%mx_mcover.OR. &
               halo_rem%i_hbeg(jjjp)+iseq_b_rem(nn)-1<1) then
             call cq_abort('index_trans: i_halo wrong ', &
                  halo_rem%i_hbeg(jjjp)+iseq_b_rem(nn)-1,gcs%mx_mcover)
          endif
          ih=halo_rem%i_halo(halo_rem%i_hbeg(jjjp)+iseq_b_rem(nn)-1)
          pairs(nr)%submat(nn)=halo_rem%ndimi(ip)*halo_rem%ndimj(ih)
          if((ip-1)*halo_rem%ni_in_halo+ih>prim%mx_iprim*halo_rem%mx_halo.OR. &
               (ip-1)*halo_rem%ni_in_halo+ih<1) then
             call cq_abort('index_trans: i_h2d wrong ', &
                  (ip-1)*halo_rem%ni_in_halo+ih,prim%mx_iprim*halo_rem%mx_halo)
          endif
          data%itran_addr(ic+nn-1)=halo_rem%i_h2d((ip-1)*halo_rem%ni_in_halo+ih)
       enddo
    enddo ! End loop over remote nodes
    return
  end subroutine index_transpose
!!***
end module trans_module
