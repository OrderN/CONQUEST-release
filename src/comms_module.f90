! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module comms_module (MPI1 version)
! ------------------------------------------------------------------------------
! Code area 2: matrices
! ------------------------------------------------------------------------------

!!****h* Conquest/comms_module *
!!  NAME
!!   comms_module
!!  PURPOSE
!!   simply holds Mquest_get and Mquest_start_send,
!!   and send_trans_data and fetch_trans_data (for the transposes) so that
!!   MPI1.2, MPI2.0 and shmem-style comms can be changed easily.
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   08/03/2002
!!  MODIFICATION HISTORY
!!   25/04/00 by DRB
!!    Changes to use generic small group types
!!   24/06/2002 dave
!!    Changed Mquest_start_send to correct (spurious) array bound errors and
!!    altered fetch_trans_data to make array copy explicit
!!   2008/02/04 17:05 dave
!!    Changed for output to file not stdout
!!   2012/01/17 17:05 dave
!!    Added routines for blip coefficient transfer (analytic blip integrals)
!1   2018/10/03 17:17 dave
!!    Changing MPI tags to conform to standard
!!  SOURCE
!!
module comms_module

  use global_module, ONLY: io_lun

  implicit none

!!***

contains

!!****f* comms_module/Mquest_start_send *
!!
!!  NAME 
!!   Mquest_start_send
!!  USAGE
!!   Mquest_start_send(comms structure,data_array,MPI requests,processor id,mx partitions on processor,sends done)
!!  PURPOSE
!!   Starts the sending of data.
!!   This subroutine is written for MPI 1.2.  It starts all the issends 
!!   that will be needed for mat_mult and get_remote.  For MPI2 it should
!!   create the windows needed, while for Cray shmems it should only 
!!   perform a barrier call 
!!  INPUTS
!! 
!! 
!!  USES
!!   mpi, datatypes, matrix_module, matrix_comms_module
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   08/03/2000
!!  MODIFICATION HISTORY
!!   24/06/2002 drb
!!    Changed the if statement (introduced get flag) to avoid array bounds checking error
!1   2018/10/03 17:17 dave
!!    Changing MPI tags to conform to standard
!!  SOURCE
!!
  subroutine Mquest_start_send(a_b_c,b,nreq,myid,mx_nponn,sends)

    ! Module usage
    use mpi
    use datatypes
    use matrix_module
    use matrix_comms_module

    implicit none

    ! Derived types
    !    type(matrix)::bmat(:)
    !    type(comms_data)::a_b_c
    type(matrix_mult)::a_b_c

    ! Passed variables
    ! Maxima
    integer :: nreq(:), myid, mx_nponn, sends
    real(double) :: b(:)
    !    integer :: bind(:)

    ! Local variables
    integer :: i,j,k,ind_part,ierr,istart,ilen1,ilen2,ilen3
    integer :: istart2,istart3,tag, lenb, lenbind, unique_parts_sent
    logical :: get

    lenb = size(b)
    lenbind = size(a_b_c%bindex)
    !write(io_lun,*) 'Lengths: ',lenb,lenbind
    sends = 0
    do i=1,a_b_c%comms%inode ! Loop over neighbour nodes
       unique_parts_sent = 0
       do j=1,a_b_c%comms%np_send(i) ! Loop over partitions to send
          ind_part = a_b_c%comms%pl_send(j,i)
          get = .false.
          if(j>1) then
             if(ind_part/=a_b_c%comms%pl_send(j-1,i)) then
                get = .true.
             end if
          else if (j==1) then
             get = .true.
          end if
          if(get) then
             unique_parts_sent = unique_parts_sent + 1
             tag = (unique_parts_sent-1)*2
             istart = a_b_c%bmat(ind_part)%array_posn
             if(istart==0) call cq_abort(' ERROR : istart zero ',ind_part)
             istart2 = a_b_c%bmat(ind_part)%nd_offset+a_b_c%bmat(ind_part)%i_nd_acc(1)
             if(istart2==0) call cq_abort(' ERROR : istart2 zero ',ind_part)
             ilen1 = a_b_c%bmat(ind_part)%n_atoms
             ilen2 = a_b_c%bmat(ind_part)%part_nabs
             ilen3 = a_b_c%bmat(ind_part)%part_nd_nabs !part_nabs
             if(istart+(3*ilen1+5*ilen2)-1>lenbind) call cq_abort('Send index error: ',istart+(3*ilen1+5*ilen2),lenbind)
             if(istart2+ilen3-1>lenb) call cq_abort('Send reals error: ',istart2+ilen3,lenb)
             ! Send the indices
             ierr = 0
             sends = sends+1
             call MPI_issend(a_b_c%bindex(istart),(3*ilen1+5*ilen2),MPI_INTEGER, &
                  a_b_c%comms%ncomm(i)-1,tag+1,MPI_COMM_WORLD,nreq(sends),ierr)
             if(ierr/=0) call cq_abort('Error sending indices ',ierr)
             ! Send xyz, sequence and elements
             if(ilen3>0)then 
                ierr = 0
                if(istart2==0) write(io_lun,*) myid,' ERROR for istart2 ',ind_part
                sends = sends+1
                call MPI_issend(b(istart2), ilen3, MPI_DOUBLE_PRECISION,&
                     a_b_c%comms%ncomm(i)-1,tag+2,MPI_COMM_WORLD,nreq(sends),ierr)
                if(ierr/=0) call cq_abort('Error sending elements ',ierr)
             endif
          endif
       enddo ! Partitions to send
    enddo ! Nodes to send to
    return
  end subroutine Mquest_start_send
!!***

! ---------------------------------------------------------------------
!   subroutine Mquest_get
! ---------------------------------------------------------------------

!!****f* comms_module/Mquest_get *
!!
!!  NAME 
!!   Mquest_get
!!  USAGE
!! 
!!  PURPOSE
!!   Calls blocking receives to get data for matrix multiplication
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
!1   2018/10/03 17:17 dave
!!    Changing MPI tags to conform to standard
!!   2018/10/04 17:31 dave
!!    Adding new argument to pass MPI tag in
!!  SOURCE
!!
  subroutine Mquest_get( mx_nponn, ilen2,ilen3,nc_part,send_node,sent_part,myid,&
       bind_rem,b_rem,lenb_rem,bind,istart,mx_babs,mx_part,tag)

    ! Module usage
    use mpi
    use datatypes
    use matrix_comms_module, ONLY: mx_msg_per_part
    use GenComms, ONLY: cq_abort

    implicit none

    ! Maxima
    integer :: mx_nponn,mx_babs,mx_part
    ! Arrays for receiving data
    integer :: bind_rem(:)
    integer :: bind(:)
    integer :: lenb_rem
    real(double) :: b_rem(lenb_rem)
    ! Miscellaneous data
    integer :: nc_part,send_node,sent_part,myid,ilen2,ilen3,size,istart,offset
    integer :: nrstat(MPI_STATUS_SIZE)
    integer :: tag

    ! Local variables
    integer :: ilen1,ierr,lenbind_rem

    !lenb_rem = size(b_rem)
    !lenbind_rem = size(bind_rem)
    ierr = 0
    ilen1 = nc_part
    !if(3*ilen1+5*ilen2>lenbind_rem) call cq_abort('Get error ',3*ilen1+5*ilen2,lenbind_rem)
    if(ilen3>lenb_rem) call cq_abort('Get error 2 ',ilen3,lenb_rem)
    call MPI_recv(bind_rem,3*ilen1+5*ilen2,MPI_INTEGER, &
         send_node-1,tag+1,MPI_COMM_WORLD,nrstat,ierr)
    if(ierr/=0) call cq_abort('Error receiving indices !',ierr)
    if(ilen3.gt.0)then ! Get xyz, sequence list and elements
       ierr = 0
       call MPI_recv(b_rem,ilen3, MPI_DOUBLE_PRECISION,send_node-1,&
            tag+2,MPI_COMM_WORLD,nrstat,ierr)
       if(ierr/=0) call cq_abort('Error receiving data !',ierr)
    endif
    return
  end subroutine Mquest_get
!!***

!!***

! ---------------------------------------------------------------------
!   subroutine Mquest_get_nonb
! ---------------------------------------------------------------------

!!****f* comms_module/Mquest_get_nonb *
!!
!!  NAME 
!!   Mquest_get_nonb
!!  USAGE
!! 
!!  PURPOSE
!!   Calls non-blocking receives to get data for matrix multiplication
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   I. Christidi
!!  CREATION DATE
!!    2023/10/06
!!  
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine Mquest_get_nonb( mx_nponn, ilen2,ilen3,nc_part,send_node,sent_part,myid,&
       bind_rem,b_rem,lenb_rem,bind,istart,mx_babs,mx_part,tag,request)

    ! Module usage
    use mpi
    use datatypes
    use matrix_comms_module, ONLY: mx_msg_per_part
    use GenComms, ONLY: cq_abort

    implicit none

    ! Maxima
    integer :: mx_nponn,mx_babs,mx_part
    ! Arrays for receiving data
    integer :: bind_rem(:)
    integer :: bind(:)
    integer :: lenb_rem
    real(double) :: b_rem(lenb_rem)
    ! Miscellaneous data
    integer :: nc_part,send_node,sent_part,myid,ilen2,ilen3,size,istart,offset
    integer :: tag
    integer :: request(2)

    ! Local variables
    integer :: ilen1,ierr,lenbind_rem

    !lenb_rem = size(b_rem)
    !lenbind_rem = size(bind_rem)
    ierr = 0
    ilen1 = nc_part
    !if(3*ilen1+5*ilen2>lenbind_rem) call cq_abort('Get error ',3*ilen1+5*ilen2,lenbind_rem)
    if(ilen3>lenb_rem) call cq_abort('Get error 2 ',ilen3,lenb_rem)
    call MPI_Irecv(bind_rem,3*ilen1+5*ilen2,MPI_INTEGER, &
         send_node-1,tag+1,MPI_COMM_WORLD,request(1),ierr)
    if(ierr/=0) call cq_abort('Error receiving indices !',ierr)
    if(ilen3.gt.0)then ! Get xyz, sequence list and elements
       ierr = 0
       call MPI_Irecv(b_rem,ilen3, MPI_DOUBLE_PRECISION,send_node-1,&
            tag+2,MPI_COMM_WORLD,request(2),ierr)
       if(ierr/=0) call cq_abort('Error receiving data !',ierr)
    endif
    return
  end subroutine Mquest_get_nonb
!!***

! ---------------------------------------------------------------------
!   subroutine send_trans_data
! ---------------------------------------------------------------------

!!****f* comms_module/send_trans_data *
!!
!!  NAME 
!!   send_trans_data
!!  USAGE
!! 
!!  PURPOSE
!!   initiate send comms needed by sbrt mat_tran.
!!   Data to be sent are matrix elements of local transpose atrans().
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
!!
!!  SOURCE
!!
  subroutine send_trans_data(mx_iprim,data,atrans,mynode, isend, nreq,pairs,size)

    use mpi
    use datatypes
    use matrix_module

    implicit none

    ! Passed variables
    integer :: mx_iprim,mynode,size

    type(trans_remote):: data
    type(pair_data) :: pairs(:)
    ! Local transpose (already done !)
    real(double) :: atrans(:)

    ! Local variables
    !ORI integer(integ) :: nnd_rem,nr,tag,n,nreq,ierr
    integer(integ) :: nnd_rem,nr,tag,n,ierr
    integer(integ),intent(out) :: nreq(:)
    integer(integ),intent(out) :: isend

    isend=0
    do nr=1,data%n_rem_node ! Loop over remote nodes
       nnd_rem=data%list_rem_node(nr)
       tag=nnd_rem-1
       n = data%n_pair(nr)
       if(nnd_rem/=mynode.AND.pairs(nr)%length>0) then
          isend=isend+1
          call MPI_issend(atrans(data%i_nd_pair_addr(nr)), pairs(nr)%length,&
               MPI_DOUBLE_PRECISION,nnd_rem-1,tag,MPI_COMM_WORLD,nreq(isend),ierr)
       endif
    enddo
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    return
  end subroutine send_trans_data
!!***

! ---------------------------------------------------------------------
!   sbrt fetch_trans_data
! ---------------------------------------------------------------------

!!****f* comms_module/fetch_trans_data *
!!
!!  NAME 
!!   fetch_trans_data
!!  USAGE
!! 
!!  PURPOSE
!!   get data elements needed by mat_tran.
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
!!   24/06/2002 drb
!!    Added bounds on the copy for a_store and atrans
!!   2008/08/12 10:46 dave
!!    Bug fix: amount copied on-processor was wrong (changed argument of source array to len not n)
!!  SOURCE
!!
  subroutine fetch_trans_data(a_store,atrans,data,len,nr,mynode,mx_iprim)

    use mpi
    use datatypes
    use matrix_module

    implicit none

    ! Passed variables
    integer :: nr,mynode,mx_iprim,len
    real(double) :: a_store(:)
    real(double) :: atrans(:)
    type(trans_remote) :: data

    ! Local variables
    integer :: tag, n, nnd_rem,ierr
    integer :: nrstat(MPI_STATUS_SIZE)

    nnd_rem=data%list_rem_node(nr)
    tag=mynode-1
    n = data%n_pair(nr)
    if(nnd_rem/=mynode.AND.len>0) then
       call MPI_recv(a_store,len,MPI_DOUBLE_PRECISION, nnd_rem-1,tag,MPI_COMM_WORLD,nrstat,ierr)
    else
       a_store(1:len) = atrans(data%i_pair_addr(nr):data%i_pair_addr(nr)+len-1)
    endif
    return
  end subroutine fetch_trans_data
!!***

! ---------------------------------------------------------------------
!   sbrt start_blip_transfer
! ---------------------------------------------------------------------

!!****f* comms_module/start_blip_transfer *
!!
!!  NAME 
!!   start_blip_transfer
!!  USAGE
!! 
!!  PURPOSE
!!   Posts non-blocking sends for blip coefficient transfers
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   2012/01/12
!!  MODIFICATION HISTORY
!!    
!!  SOURCE
!!
  subroutine start_blip_transfer(nreq,sends,mxng)

    use matrix_data, ONLY: blip_trans
    use mpi
    use GenComms, ONLY: cq_abort, myid
    use support_spec_format, ONLY : coefficient_array
    
    implicit none

    ! Passed
    integer :: nreq(:)
    integer :: sends,mxng

    ! Local
    integer :: i,j,tag, ind_part, ierr
    logical :: get

    sends = 0
    nreq = 0
    do i=1,blip_trans%nproc  ! Loop over neighbour nodes
       do j=1,blip_trans%np_send(i) ! Loop over partitions to send
          ind_part = blip_trans%pl_send(j,i)
          get = .false.
          if(j>1) then
             if(ind_part/=blip_trans%pl_send(j-1,i)) then
                get = .true.
             end if
          else if (j==1) then
             get = .true.
          end if
          if(get.AND.blip_trans%partlen(ind_part)>0) then
             sends = sends + 1
             ! This defines a unique tag for each message sent by this node
             tag = (blip_trans%ncomm(i)-1)*mxng + ind_part
             call MPI_issend(coefficient_array(blip_trans%partst(ind_part)),blip_trans%partlen(ind_part),MPI_DOUBLE_PRECISION, &
                  blip_trans%ncomm(i)-1,tag,MPI_COMM_WORLD,nreq(sends),ierr)
             if(ierr/=0) call cq_abort('Error sending blip coefficients for part: ',j,ierr)
          end if
       end do
    end do
    return
  end subroutine start_blip_transfer
!!***

! ---------------------------------------------------------------------
!   sbrt fetch_blips
! ---------------------------------------------------------------------

!!****f* comms_module/fetch_blips *
!!
!!  NAME 
!!   fetch_blips
!!  USAGE
!! 
!!  PURPOSE
!!   Fetches blip coefficients for a partition
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   2012/01/12
!!  MODIFICATION HISTORY
!!    
!!  SOURCE
!!
  subroutine fetch_blips(part_blips,part_len,send_proc,tag)
  
    use mpi
    use datatypes
    use GenComms, ONLY: cq_abort, myid
    use support_spec_format, ONLY : coefficient_array

    implicit none
  
    ! Passed
    integer :: part_len, send_proc, tag
    real(double), dimension(part_len) :: part_blips
  
    ! Local
    integer :: ierr, stat
    integer :: nrstat(MPI_STATUS_SIZE)
  
    ierr = 0
    !tag = myid
    if(part_len>0) then!send_proc/=myid+1.AND.part_len>0) then
       call MPI_recv(part_blips,part_len,MPI_DOUBLE_PRECISION,send_proc, tag, MPI_COMM_WORLD, nrstat, ierr) 
       if(ierr/=0) call cq_abort('Error receiving blip coefficients: ',ierr)
    end if
    return
  end subroutine fetch_blips
!!***

end module comms_module
