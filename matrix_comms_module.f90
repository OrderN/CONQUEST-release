! $Id$
! -----------------------------------------------------------
! Module matrix_comms_module
! -----------------------------------------------------------
! Code area 2: Matrices
! -----------------------------------------------------------

!!****h* Conquest/matrix_comms_module
!!  NAME
!!   matrix_comms_module
!!  PURPOSE
!!   module matrix_comms_module: contains all the bits of communication 
!!   required for matrix comms: Mquest_get and start_send, along with 
!!   end_part_comms and init_part_comms
!! 
!!   Actually, the comms parts are now in comms_module which can be
!!   swapped in and out depending on what communications protocol is
!!   used - MPI1, MPI2, shmem etc.
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   25/11/99
!!  MODIFICATION HISTORY
!!   7/12/99 by D.R.Bowler
!!    Memory allocated to part_list in init_mult_comms corrected
!!   20/06/2001 dave
!!    Added ROBODoc header, RCS Id and Log tags and removed stops
!!    Added GenComms to init_mult_comms
!!   2008/02/06 08:22 dave
!!    Changed for output to file not stdout
!!  SOURCE
!!
module matrix_comms_module

  ! Module usage
  use global_module

  implicit none

  integer, parameter :: mx_msg_per_part=10 ! max msgs sent/recv per part
!!***

contains

!!****f* matrix_comms_module/end_part_comms *
!!
!!  NAME 
!!   end_part_comms
!!  USAGE
!! 
!!  PURPOSE
!!   Builds ibpart_rem from the offsets contained in npxyz_rem
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   25/11/99 ?
!!  MODIFICATION HISTORY
!!   20/06/2001 dave
!!    Added ROBODoc header
!!  SOURCE
!!
  subroutine end_part_comms(myid,n_cont, &
       nbnab_rem,ibind_rem,npxyz_rem,ibpart_rem, &
       ncover_yz,ncoverz)

    ! Module usage
    use datatypes
    use GenComms, ONLY: cq_abort

    ! Passed variables
    integer :: myid
    integer :: n_cont,ncover_yz,ncoverz
    integer(integ) :: nbnab_rem(:)
    integer(integ) :: ibind_rem(:)
    integer(integ) :: ibpart_rem(:)
    integer(integ) :: npxyz_rem(:)

    ! Local variables
    integer :: noff,ierr,ni,k, len1, len2

    len1 = size(ibpart_rem)
    len2 = size(npxyz_rem)
    ! Now build the partition list from the offsets
    noff=0
    ! Loop over atoms in partition and their neighbours
    do k=1,n_cont
       do ni=1,nbnab_rem(k)
          if(ibind_rem(k)+ni-1>len1) call cq_abort('Overflow in end_part_comms part: ',ibind_rem(k)+ni-1,len1)
          if(noff+3+3*(ni-1)>len2) call cq_abort('Overflow in end_part_comms part: ',noff+3+3*(ni-1),len1)
          ibpart_rem(ibind_rem(k)+ni-1)= &
               npxyz_rem(noff+1+3*(ni-1))*ncover_yz &
               +npxyz_rem(noff+2+3*(ni-1))*ncoverz &
               +npxyz_rem(noff+3+3*(ni-1))
       enddo
       noff = noff+3*nbnab_rem(k)
    enddo
    return
  end subroutine end_part_comms
!!***

!!****f* matrix_comms_module/Mquest_get_local *
!!
!!  NAME 
!!   Mquest_get_local
!!  USAGE
!! 
!!  PURPOSE
!!   Performs a copy when the B data for multiplication is local
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   25/11/99 ?
!!  MODIFICATION HISTORY
!!   20/06/2001 dave
!!    Added ROBODoc header
!!  TODO
!!   Change ndims to maxima so that we can use dcopy 15/06/2001 dave
!!  SOURCE
!!
  subroutine Mquest_get_local(ip, &
       bind_rem,b_rem,lenb_rem,bind,bmat,&
       ind_part,b,myid)

    use datatypes
    use matrix_module

    implicit none

    integer :: ip,ilen3,istart,istart2,istart3,ii,ind_part,jj,kk
    integer :: ni,myid,indsize

    integer :: offset, lenb_rem
    integer(integ), dimension(:) :: bind_rem
    integer(integ), dimension(:) :: bind

    type(matrix), dimension(:) :: bmat

    real(double) :: b_rem(lenb_rem)
    real(double) :: b(:)
    integer :: lenb, lenbind,lenbind_rem

    !lenb = size(b)
    !lenb_rem = size(b_rem)
    !lenbind = size(bind)
    !lenbind_rem = size(bind_rem)
    ! Attributes for index copy
    indsize = (3*bmat(ip)%n_atoms+5*bmat(ip)%part_nabs)
    offset = bmat(ip)%array_posn-1
    !if(indsize>lenbind_rem) call cq_abort('bind_rem overrun: ',indsize,lenbind_rem)
    !if(offset+indsize>lenbind) call cq_abort('bind overrun: ',offset+indsize,lenbind)
    bind_rem(1:indsize) = bind(offset+1:offset+indsize)
    ! Attributes for elements copy
    ilen3 = bmat(ip)%part_nd_nabs
    istart2 = bmat(ip)%nd_offset+bmat(ip)%i_nd_acc(1)-1
    if(ilen3>lenb_rem) call cq_abort('b_rem overrun: ',ilen3,lenb_rem)
    !if(istart2+ilen3>lenb) call cq_abort('b overrun: ',ilen3+istart2,lenb)
    b_rem(1:ilen3) = b(istart2+1:istart2+ilen3)
    return
  end subroutine Mquest_get_local
!!***

!!****f* matrix_comms_module/init_mult_comms *
!!
!!  NAME 
!!   init_mult_comms
!!  USAGE
!! 
!!  PURPOSE
!!   sbrt init_mult_comms: initiates the multiplication between the matrices 
!!   given by ahalo and bmat.  It calculates all the communications 
!!   which need to be done, and stores all the data in the derived type a_b_c
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   25/11/99 
!!  MODIFICATION HISTORY
!!   20/06/2001 dave
!!    Added ROBODoc header
!!    Added GenComms
!!   2008/02/11 16:55 dave
!!    Added MPI_Wait calls for non-blocking sends
!!   2011/12/06 17:05 dave
!!    Bug fix for arrays passed as nrstat
!!  SOURCE
!!
  subroutine init_mult_comms(parts,a_b_c,myid)

    ! Module usage
    use mpi
    use datatypes
    use global_module
    use matrix_module
    use basic_types
    use GenComms, ONLY: cq_abort
    use maxima_module, ONLY: maxnabaprocs

    implicit none

    ! Passed variables
    type(matrix_mult)  ::a_b_c
    type(group_set)  :: parts

    integer :: ilen2(parts%mx_ngonn,numprocs)
    integer :: ilen3(parts%mx_ngonn,numprocs)
    integer :: istart(parts%mx_ngonn,numprocs)
    integer :: myid

    ! Local variables
    integer :: ierr,kpart,ind_part,i,j,nod,stat
    integer :: nreq(3*maxnabaprocs)
    integer :: nrstat(MPI_STATUS_SIZE,3*maxnabaprocs)
    integer :: np_get(maxnabaprocs)  
    integer, pointer :: part_list(:,:)

    logical :: new

    ! Allocate memory
    call allocate_comms_data(a_b_c%comms,a_b_c%ahalo%mx_part, &
         maxnabaprocs,parts%mx_ngonn)
    allocate(part_list(a_b_c%ahalo%mx_part,maxnabaprocs),STAT=stat)
    if(stat/=0) then
       call cq_abort('init_mult_comms: error alloc part_list')
    endif
    ! Make a list of which nodes I need to go to to get data (and hence which
    ! will be coming to me)
    np_get = 0
    part_list = 0
    ilen2 = 0
    ilen3 = 0
    istart = 0
    a_b_c%comms%ncomm = 0
    a_b_c%comms%inode = 0
    ! Now loop over partitions and work out where to get data from/send data to
    if(a_b_c%ahalo%np_in_halo>a_b_c%ahalo%mx_part) call cq_abort('ERROR in HALO',a_b_c%ahalo%np_in_halo,a_b_c%ahalo%mx_part)
    do kpart=1,a_b_c%ahalo%np_in_halo
       new = .true.
       ind_part=a_b_c%ahalo%lab_hcell(kpart)
       nod = parts%i_cc2node(ind_part)
       if(a_b_c%comms%inode > 0) then
          do i=1,a_b_c%comms%inode  ! is this unnecessary (kpart is in node order)
             if(nod.eq.a_b_c%comms%ncomm(i)) then  !If we've seen this node before
                new = .false.         
                np_get(i) = np_get(i)+1
                part_list(np_get(i),i) = parts%i_cc2seq(ind_part)
                a_b_c%comms%neigh_node_list(kpart) = i
             endif
          enddo
       end if ! (a_b_c%comms%inode > 0)
       if(new.AND.nod.ne.myid+1)then 
          a_b_c%comms%inode = a_b_c%comms%inode + 1
          a_b_c%comms%ncomm(a_b_c%comms%inode) = nod
          np_get(a_b_c%comms%inode) = np_get(a_b_c%comms%inode)+1
          part_list(np_get(a_b_c%comms%inode),a_b_c%comms%inode) &
               = parts%i_cc2seq(ind_part)
          a_b_c%comms%neigh_node_list(kpart) = a_b_c%comms%inode
       else if(new.AND.nod.eq.myid+1)then
          a_b_c%comms%neigh_node_list(kpart) = a_b_c%comms%inode
       endif
    enddo    
    ! Now send to and receive from those nodes our comms lists
    ! We post non-blocking sends for all our data, and then get
    ! what we need.  This should be guaranteed to work under MPI.
    if(a_b_c%comms%inode > 0) then
       do i=1,a_b_c%comms%inode
          call MPI_issend(np_get(i),1,MPI_INTEGER,a_b_c%comms%ncomm(i)-1,&
               1,MPI_COMM_WORLD,nreq(i*2-1),ierr)
          call MPI_issend(part_list(1,i),np_get(i),&
               MPI_INTEGER,a_b_c%comms%ncomm(i)-1,2,MPI_COMM_WORLD,nreq(i*2),ierr)
       enddo
       ! Now do the blocking gets - we need to have this data
       do i=1,a_b_c%comms%inode
          call MPI_recv(a_b_c%comms%np_send(i),1,MPI_INTEGER, &
               a_b_c%comms%ncomm(i)-1,1,MPI_COMM_WORLD,nrstat(1:MPI_STATUS_SIZE,i*2-1),ierr)
          call MPI_recv(a_b_c%comms%pl_send(1,i),a_b_c%comms%np_send(i),&
               MPI_INTEGER,a_b_c%comms%ncomm(i)-1,2,&
               MPI_COMM_WORLD,nrstat(1:MPI_STATUS_SIZE,i*2),ierr)
       enddo
       ! Add MPI_Wait calls
       do i=1,a_b_c%comms%inode
          call MPI_Wait(nreq(i*2-1),nrstat(1:MPI_STATUS_SIZE,1),ierr)
          call MPI_Wait(nreq(i*2),nrstat(1:MPI_STATUS_SIZE,2),ierr)
       end do
       ! The next two send/recv pairs are for the ilen2 parameter - this
       ! is the sum of the number of B neighbours for all atoms in a 
       ! partition.  It's rather useful.
       ! -------------------------------------------------------
       ! First, work out what we need to send
       do i=1,a_b_c%comms%inode
          do kpart=1,a_b_c%comms%np_send(i)
             if(kpart.gt.1) then
                if(a_b_c%comms%pl_send(kpart,i)/=a_b_c%comms%pl_send(kpart-1,i)) then
                   ilen2(a_b_c%comms%pl_send(kpart,i),i) = &
                        a_b_c%bmat(a_b_c%comms%pl_send(kpart,i))%part_nabs
                   ilen3(a_b_c%comms%pl_send(kpart,i),i) = &
                        a_b_c%bmat(a_b_c%comms%pl_send(kpart,i))%part_nd_nabs
                   istart(a_b_c%comms%pl_send(kpart,i),i) = &
                        a_b_c%bmat(a_b_c%comms%pl_send(kpart,i))%nd_offset+ &
                        a_b_c%bmat(a_b_c%comms%pl_send(kpart,i))%i_nd_acc(1)
                endif
             else
                ilen2(a_b_c%comms%pl_send(kpart,i),i) = &
                     a_b_c%bmat(a_b_c%comms%pl_send(kpart,i))%part_nabs
                ilen3(a_b_c%comms%pl_send(kpart,i),i) = &
                     a_b_c%bmat(a_b_c%comms%pl_send(kpart,i))%part_nd_nabs
                istart(a_b_c%comms%pl_send(kpart,i),i) = &
                     a_b_c%bmat(a_b_c%comms%pl_send(kpart,i))%nd_offset+ &
                     a_b_c%bmat(a_b_c%comms%pl_send(kpart,i))%i_nd_acc(1)
             endif
          enddo
       enddo
       ! Now send our stuff, and then get what we need
       do i=1,a_b_c%comms%inode
          call MPI_issend(ilen2(1,i),parts%mx_ngonn,MPI_INTEGER, &
               a_b_c%comms%ncomm(i)-1,3,MPI_COMM_WORLD,nreq(i*3-2),ierr)
          call MPI_issend(ilen3(1,i),parts%mx_ngonn,MPI_INTEGER, &
               a_b_c%comms%ncomm(i)-1,4,MPI_COMM_WORLD,nreq(i*3-1),ierr)
          call MPI_issend(istart(1,i),parts%mx_ngonn,MPI_INTEGER, &
               a_b_c%comms%ncomm(i)-1,5,MPI_COMM_WORLD,nreq(i*3),ierr)
       enddo
       do i=1,a_b_c%comms%inode
          call MPI_recv(a_b_c%comms%ilen2rec(1,i),parts%mx_ngonn,&
               MPI_INTEGER,a_b_c%comms%ncomm(i)-1,3,&
               MPI_COMM_WORLD,nrstat(1:MPI_STATUS_SIZE,i*2),ierr)
          call MPI_recv(a_b_c%comms%ilen3rec(1,i),parts%mx_ngonn,&
               MPI_INTEGER,a_b_c%comms%ncomm(i)-1,4,&
               MPI_COMM_WORLD,nrstat(1:MPI_STATUS_SIZE,i*2),ierr)
          call MPI_recv(a_b_c%comms%istart(1,i),parts%mx_ngonn,&
               MPI_INTEGER,a_b_c%comms%ncomm(i)-1,5,&
               MPI_COMM_WORLD,nrstat(1:MPI_STATUS_SIZE,i*2),ierr)
       enddo
       do i=1,a_b_c%comms%inode
          call MPI_Wait(nreq(i*3-2),nrstat(1:MPI_STATUS_SIZE,1),ierr)
          call MPI_Wait(nreq(i*3-1),nrstat(1:MPI_STATUS_SIZE,2),ierr)
          call MPI_Wait(nreq(i*3),nrstat(1:MPI_STATUS_SIZE,3),ierr)
       end do
    end if !(a_b_c%comms%inode > 0)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    deallocate(part_list,STAT=stat)
    if(stat/=0) then
       call cq_abort('init_mult_comms: error dealloc part_list')
    endif
    return
  end subroutine init_mult_comms
!!***

!!****f* matrix_comms_module/find_neighbour_procs *
!!
!!  NAME 
!!   find_neighbour_nodes
!!  USAGE
!! 
!!  PURPOSE
!!   Finds maximum nodes we will exchange multiply data with, based on halo of LSL
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   2006/09/07
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine find_neighbour_procs(parts,halo)

    use basic_types, ONLY: group_set
    use matrix_module, ONLY: matrix_halo
    use maxima_module, ONLY: maxnabaprocs
    use GenComms, ONLY: myid, cq_abort
    use global_module, ONLY: numprocs

    implicit none

    ! Passed variables
    type(matrix_halo) :: halo
    type(group_set)  :: parts

    ! Local variables
    integer :: kpart, ind_part, nod,i
    integer, dimension(numprocs) :: listnodes
    logical :: new

    listnodes = 0
    maxnabaprocs = 0
    do kpart=1,halo%np_in_halo
       new = .true.
       ind_part=halo%lab_hcell(kpart)
       nod = parts%i_cc2node(ind_part)
       if(maxnabaprocs>0) then
          do i=1,maxnabaprocs
             if(nod==listnodes(i)) new = .false.         
          enddo
       end if !(maxnabaprocs>0)
       if(new.AND.nod.ne.myid+1)then 
          maxnabaprocs = maxnabaprocs + 1
          listnodes(maxnabaprocs) = nod
       endif
    enddo
    if(maxnabaprocs==0) maxnabaprocs = 1
    if(iprint_mat>2.AND.myid==0) write(io_lun,fmt='(2x,"Proc: ",i5," has ",i5," neighbour processors")') myid+1,maxnabaprocs
  end subroutine find_neighbour_procs
!!***
end module matrix_comms_module
