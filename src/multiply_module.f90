! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module multiply_module
! ------------------------------------------------------------------------------
! Code area 2: matrices
! ------------------------------------------------------------------------------

!!****h* Conquest/multiply_module *
!!  NAME
!!   multiply_module
!!  PURPOSE
!!   Contains subroutines related to matrix manipulation:
!!   multiplication (plus local transpose, prefetch, kernel etc);
!!   transposition; adding and copying different range matrices.
!!  USES
!!   GenComms, basic_types, comms_module, datatypes, global_module, matrix_comms_module,
!!   matrix_module, mpi, trans_module
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   25/11/99
!!  MODIFICATION HISTORY
!!   7/12/99 by D.R.Bowler
!!    Sizes of the _rem variables in mat_mult changed from mx_part*mx_bnab
!!    to mx_part*mx_babs (subtle difference which was causing array overflow)
!!    Also changed parameters passed for these sizes to end_part_comms etc
!!   27/01/00 by D.R.Bowler
!!    Now uses new comms: all indices in the matrix type are stored by
!!    partition, and point to different areas of an array; that array is
!!    all that is sent when we transfer the indices
!!   20/06/2001 dave
!!    Added ROBODoc headers, RCS Id and Log tags and cq_abort throughout
!!   24/04/2002 dave
!!    Corrected ROBODoc headers and added RCS static Id tag, and USES field
!!   17/06/2002 dave
!!    Tiny change to headers
!!   20/06/2002 drb
!!    Added numbers to cq_abort calls in loc_trans to make them more meaningful
!!   20/06/2002 dave
!!    Changed dimension of nreqs (away from 500 to actual maximum) and made atrans dynamically allocated rather
!!    than automatically allocated.
!!   10:59, 13/02/2006 drb
!!    Removed prune, graft routines (no longer used) and small tidying in preparation for variable NSF
!!   2008/02/06 08:26 dave
!!    Changed for output to file not stdout
!!   2008/05/22 ast
!!    Added timers
!!   2012/08/20 L.Tong
!!   - Moved multiplication kernels to separate modules
!!     multiply_kernel_VERSION.f90
!!   2014/09/15 18:30 lat
!!    fixed call start/stop_timer to timer_module (not timer_stdlocks_module !)
!!  SOURCE
!!
module multiply_module

  use global_module,          only: io_lun
  use timer_module,           only: start_timer, stop_timer
  use timer_stdclocks_module, only: tmr_std_allocation

  implicit none

!!***

contains

  !!****f* multiply_module/mat_mult *
  !!
  !!  NAME
  !!   mat_mult - sparse matrix multiplication
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Multiply together a and b and produce c.
  !!   There are two types of multiplication - "extension" where
  !!   Rc is close to Ra + Rb, and "reduction" where Rc is much
  !!   smaller.  In practice, the latter case is actually done as
  !!   A = C.B.
  !!
  !!   See the notes "A comprehensive scheme for matrix multiplication".
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
  !!    Added ROBODoc header and cq_abort and tidied
  !!   20/06/2002 dave
  !!    Changed dimension of nreqs to be mx_nponn*mx_neighbour_nodes
  !!    (the actual maximum that will be required) and converted
  !!    atrans to being dynamically allocated rather than automatic.
  !!   10:44, 06/03/2003 drb
  !!    Changed size of nreqs (temporary fix)
  !!   09:19, 11/05/2005 dave
  !!    Removed unnecessary MPI_Wait (second one on nreqs)
  !!   2008/05/22 ast
  !!    Added timers
  !!   2012/08/20 L.Tong
  !!   - Removed conditionals for choosing GEMM or normal kernels.
  !!     The choice will be done during compilation time.
  !!   2012/09/05 L.Tomg
  !!   - Added timer for matrix multiplication
  !!   2018/10/04 17:26 dave
  !!    Added counter for partitions received from different processes (nodes)
  !!    to use for tags to be compliant with MPI standard
  !!  SOURCE
  !!
  subroutine mat_mult(myid,a,lena,b,lenb,c,lenc,a_b_c,debug)

    ! Module usage
    use datatypes
    use numbers
    use global_module
    use matrix_module
    use basic_types
    use matrix_comms_module
    use comms_module
    use GenComms, only: my_barrier, cq_abort, mtime
    use mpi
    use multiply_kernel
    use timer_stdclocks_module, only: tmr_std_matmult

    implicit none

    ! Mult type - contains indices, small groups etc
    type(matrix_mult)::a_b_c
    ! Matrices
    integer :: lena, lenb, lenc
    real(double) :: a(lena)
    real(double) :: b(lenb)
    real(double) :: c(lenc)
    ! Node id - should be in a module ?
    integer :: myid
    integer, OPTIONAL :: debug

    ! Local variables
    ! This will be dynamically allocated/deallocated by the system
    real(double), allocatable, dimension(:) :: atrans
    integer :: lab_const
    integer :: invdir,ierr,kpart,ind_part,ncover_yz,n_which,ipart,nnode
    integer :: icall,n_cont,kpart_next,ind_partN,k_off
    integer :: stat,ilen2,lenb_rem
    ! Remote variables to be allocated
    integer(integ),allocatable :: ibpart_rem(:)
    real(double),allocatable :: b_rem(:)
    ! Remote variables which will point to part_array
    integer(integ),pointer :: nbnab_rem(:)
    integer(integ),pointer :: ibseq_rem(:)
    integer(integ),pointer :: ibind_rem(:)
    integer(integ),pointer :: ib_nd_acc_rem(:)
    integer(integ),pointer :: npxyz_rem(:)
    integer(integ),pointer :: ibndimj_rem(:)
    ! Arrays for remote variables to point to
    integer, target :: part_array(3*a_b_c%parts%mx_mem_grp+ &
         5*a_b_c%parts%mx_mem_grp*a_b_c%bmat(1)%mx_abs)
    integer, dimension(:), allocatable :: nreqs
    integer :: offset,sends,i,j
    integer, dimension(MPI_STATUS_SIZE) :: mpi_stat
    integer, allocatable, dimension(:) :: recv_part
    real(double) :: t0,t1

    logical :: new_partition

    call start_timer(tmr_std_matmult)

    call start_timer(tmr_std_allocation)
    if(iprint_mat>3.AND.myid==0) t0 = mtime()
    ! Allocate memory for the elements
    allocate(ibpart_rem(a_b_c%parts%mx_mem_grp*a_b_c%bmat(1)%mx_abs),STAT=stat)
    if(stat/=0) call cq_abort('mat_mult: error allocating ibpart_rem')
    !allocate(atrans(a_b_c%amat(1)%length),STAT=stat)
    allocate(atrans(lena),STAT=stat)
    if(stat/=0) call cq_abort('mat_mult: error allocating atrans')
    allocate(recv_part(0:a_b_c%comms%inode),STAT=stat)
    if(stat/=0) call cq_abort('mat_mult: error allocating recv_part')
    recv_part = zero
    call stop_timer(tmr_std_allocation)
    !write(io_lun,*) 'Sizes: ',a_b_c%comms%mx_dim3*a_b_c%comms%mx_dim2*a_b_c%parts%mx_mem_grp*a_b_c%bmat(1)%mx_abs,&
    !     a_b_c%parts%mx_mem_grp*a_b_c%bmat(1)%mx_abs,a_b_c%comms%mx_dim3*a_b_c%comms%mx_dim1* &
    !     a_b_c%prim%mx_iprim*a_b_c%amat(1)%mx_nab
    sends = 0
    do i=1,a_b_c%comms%inode
       if(a_b_c%comms%np_send(i)>0) then
          do j=1,a_b_c%comms%np_send(i)
             sends = sends+1
          end do
       end if
    end do
    call start_timer(tmr_std_allocation)
    allocate(nreqs(sends*2),STAT=stat)
    if(stat/=0) call cq_abort("mat_mult: Error allocating nreqs",sends,stat)
    call stop_timer(tmr_std_allocation)
    sends = 0
    ! For type 1, form local transpose of my data (a_halo) and zero C-matrix
    invdir=0
    !write(io_lun,*) 'Local trans'
    if(a_b_c%mult_type.eq.1) then
       call loc_trans( a_b_c%ltrans, a_b_c%ahalo,a,lena,atrans,lena,invdir)
       c = zero
    else if(a_b_c%mult_type.eq.2) then
       atrans = zero
    end if
    ! Start the send procedure
    !write(io_lun,*) 'Sending'
    call Mquest_start_send(a_b_c,b,nreqs,myid,a_b_c%prim%mx_ngonn,sends)
    !write(io_lun,*) 'Returned ',a_b_c%ahalo%np_in_halo,myid
    ncover_yz=a_b_c%gcs%ncovery*a_b_c%gcs%ncoverz

    !$omp parallel default(shared)
    main_loop: do kpart = 1,a_b_c%ahalo%np_in_halo

       !$omp master
       icall=1
       ind_part = a_b_c%ahalo%lab_hcell(kpart)
       new_partition = .true.
       
       ! Check if this is a periodic image of the previous partition
       if(kpart>1) then
          if(ind_part.eq.a_b_c%ahalo%lab_hcell(kpart-1)) then
             new_partition = .false.
          end if
       end if

       if(new_partition) then
          ! Get the data
          ipart = a_b_c%parts%i_cc2seq(ind_part)
          nnode = a_b_c%comms%neigh_node_list(kpart)
          recv_part(nnode) = recv_part(nnode)+1
          if(allocated(b_rem)) deallocate(b_rem)
          if(a_b_c%parts%i_cc2node(ind_part)==myid+1) then
             lenb_rem = a_b_c%bmat(ipart)%part_nd_nabs
          else
             lenb_rem = a_b_c%comms%ilen3rec(ipart,nnode)
          end if
          allocate(b_rem(lenb_rem))
          call prefetch(kpart,a_b_c%ahalo,a_b_c%comms,a_b_c%bmat,icall,&
               n_cont,part_array,a_b_c%bindex,b_rem,lenb_rem,b,myid,ilen2,&
               mx_msg_per_part,a_b_c%parts,a_b_c%prim,a_b_c%gcs,(recv_part(nnode)-1)*2)
          ! Now point the _rem variables at the appropriate parts of
          ! the array where we will receive the data
          offset = 0
          nbnab_rem => part_array(offset+1:offset+n_cont)
          offset = offset+n_cont
          ibind_rem => part_array(offset+1:offset+n_cont)
          offset = offset+n_cont
          ib_nd_acc_rem => part_array(offset+1:offset+n_cont)
          offset = offset+n_cont
          ibseq_rem => part_array(offset+1:offset+ilen2)
          offset = offset+ilen2
          npxyz_rem => part_array(offset+1:offset+3*ilen2)
          offset = offset+3*ilen2
          ibndimj_rem => part_array(offset+1:offset+ilen2)
          if(offset+ilen2>3*a_b_c%parts%mx_mem_grp+ &
               5*a_b_c%parts%mx_mem_grp*a_b_c%bmat(1)%mx_abs) then
             call cq_abort('mat_mult: error pointing to part_array ',kpart)
          end if
          ! Create ibpart_rem
          call end_part_comms(myid,n_cont,nbnab_rem,ibind_rem,npxyz_rem,&
               ibpart_rem,ncover_yz,a_b_c%gcs%ncoverz)
       end if
       
       k_off=a_b_c%ahalo%lab_hcover(kpart) ! --- offset for pbcs
       ! Omp master doesn't include a implicit barrier. We want master
       ! to be finished with comms before calling the multiply kernels
       ! hence the explicit barrier
       !$omp end master
       !$omp barrier
       
       if(a_b_c%mult_type.eq.1) then  ! C is full mult
          call m_kern_max( k_off,kpart,ib_nd_acc_rem, ibind_rem,nbnab_rem,&
               ibpart_rem,ibseq_rem,ibndimj_rem,&
               atrans,b_rem,c,a_b_c%ahalo,a_b_c%chalo,a_b_c%ltrans,&
               a_b_c%bmat(1)%mx_abs,a_b_c%parts%mx_mem_grp, &
               a_b_c%prim%mx_iprim, lena, lenb_rem, lenc)
       else if(a_b_c%mult_type.eq.2) then ! A is partial mult
          call m_kern_min( k_off,kpart,ib_nd_acc_rem, ibind_rem,nbnab_rem,&
               ibpart_rem,ibseq_rem,ibndimj_rem,&
               atrans,b_rem,c,a_b_c%ahalo,a_b_c%chalo,a_b_c%ltrans,&
               a_b_c%bmat(1)%mx_abs,a_b_c%parts%mx_mem_grp, &
               a_b_c%prim%mx_iprim, lena, lenb_rem, lenc)
       end if
       !$omp barrier
    end do main_loop
    !$omp end parallel
    call start_timer(tmr_std_allocation)
    if(allocated(b_rem)) deallocate(b_rem)
    call stop_timer(tmr_std_allocation)
    ! --------------------------------------------------
    ! End of the main loop over partitions in the A-halo
    ! --------------------------------------------------
    ! sends is only incremented in the MPI version of Mquest_start_send, so this works
    !write(io_lun,*) 'Done loop ',myid,sends
    if(sends>0) then
       do i=1,sends
          call MPI_Wait(nreqs(i),mpi_stat,ierr)
          if(ierr/=0) call cq_abort("Error waiting for send to finish",i)
          !write(io_lun,*) 'Send done ',i,myid
       end do
    end if
    call my_barrier
    call start_timer(tmr_std_allocation)
    deallocate(nreqs,STAT=stat)
    if(stat/=0) call cq_abort('mat_mult: error deallocating nreqs')
    call stop_timer(tmr_std_allocation)
    call my_barrier
    ! --- for type 2, make backward local transpose of A-matrix -----------
    if(a_b_c%mult_type.eq.2) then
       invdir=1
       call loc_trans( a_b_c%ltrans, a_b_c%ahalo,a,lena,atrans,lena,invdir)
    end if
    call my_barrier
    call start_timer(tmr_std_allocation)
    deallocate(atrans,STAT=stat)
    if(stat/=0) call cq_abort('mat_mult: error deallocating atrans')
    call stop_timer(tmr_std_allocation)
    call my_barrier
    call start_timer(tmr_std_allocation)
    deallocate(ibpart_rem,STAT=stat)
    if(stat/=0) call cq_abort('mat_mult: error deallocating ibpart_rem')
    deallocate(recv_part,STAT=stat)
    if(stat/=0) call cq_abort('mat_mult: error deallocating recv_part')
    call stop_timer(tmr_std_allocation)
    call my_barrier
    !deallocate(b_rem,STAT=stat)
    !if(stat/=0) call cq_abort('mat_mult: error deallocating b_rem')
    !call my_barrier
    if(iprint_mat>3.AND.myid==0) then
       t1 = mtime()
       write(io_lun,fmt='(10x,a,f16.6)') 'mult time: ',t1-t0
    end if

    call stop_timer(tmr_std_matmult)

    return
  end subroutine mat_mult
  !!***


  !!****f* multiply_module/mat_trans *
  !!
  !!  NAME
  !!   mat_trans - global transpose of a matrix
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Takes the global transpose of a matrix.  Many of the preparatory
  !!   and subsidiary routines for this are in trans_module.
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
  !!    Added ROBODoc header and cq_abort
  !!  SOURCE
  !!
  subroutine mat_trans(nnd,at_rem,at,ahalo,a,atrans,ltrans,mx_iprim)

    use datatypes
    use matrix_module
    use comms_module
    !MPI_Wait  T. Miyazaki 2003/11/18
    !use mpi
    use GenComms, ONLY: cq_abort

    implicit none

    integer(integ)     :: nnd,mx_iprim,dr
    type(trans_remote) :: at_rem
    type(matrix_trans) :: at
    type(matrix_halo)  :: ahalo
    real(double)       :: a(:),atrans(:)

    real(double)       :: ltrans(:)

    !MPI_Wait  T. Miyazaki 2003/11/18
    !integer  :: nreqs(1:1000), isend, ii, ierr=0
    !integer, dimension(MPI_STATUS_SIZE) :: mpi_stat

    ! Form the local transpose of A
    atrans = 0.0_double
    ltrans = 0.0_double
    !call loc_trans(at,ahalo,a,ltrans,0)
    ! Send out the data
    !call send_trans_data(mx_iprim,at_rem,ltrans,nnd, isend, nreqs)
    ! Receive it and built AT
    !call mat_tran(mx_iprim,at%mx_nab,nnd,at_rem,atrans,ltrans)

  !MPI_Wait  T. Miyazaki 2003/11/18
!    do ii= 1, isend
!      call MPI_Wait(nreqs(ii),mpi_stat,ierr)
!      if(ierr/=0) call cq_abort("Error in mat_trans: waiting for send to finish",ii)
!    end do
   return
  end subroutine mat_trans
  !!***


  !!****f* multiply_module/loc_trans *
  !!
  !!  NAME
  !!   loc_trans - local (to processor) transpose
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Performs a local transpose, so that A_ik (indexed as primary set
  !!   atoms, halo atoms) becomes A_ki (halo atoms, primary set atoms).
  !!   This is used extensively by the matrix multiplication.
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
  !!    Added ROBODoc header and cq_abort
  !!   20/06/2002 drb
  !!    Added numbers to cq_abort calls so that they are more meaningful
  !!  SOURCE
  !!
  subroutine loc_trans(at,ahalo,a,lena,atrans,lenat,invdir)

    use datatypes
    use GenComms, ONLY: cq_abort
    use maxima_module, ONLY: maxatomsproc
    use matrix_module, ONLY: matrix_trans, matrix_halo
    use basic_types, ONLY: primary_set
    use primary_module, ONLY: bundle

    implicit none

    type(matrix_trans) :: at
    type(matrix_halo)  :: ahalo

    integer :: invdir,irc,ierr,ni,nn
    integer :: i,istart,n1,n2,tpos,nd1,nd2
    ! ---------------------------------------------------------------------
    integer :: lena,lenat
    real(double) :: a(lena)
    real(double) :: atrans(lenat)

    ! --- invdir=0: normal -> retro; invdir=1: retro -> normal ------------
    if(invdir.eq.0) then
       do ni=1,ahalo%ni_in_halo
          !write(io_lun,*) 'Atom no: ',ni
          if(at%n_hnab(ni)<=0) then
             call cq_abort('loc_trans: no. of neighbours < 1')
          end if
          if(at%i_beg(ni)+at%n_hnab(ni)-1>maxatomsproc*at%mx_nab) then
             call cq_abort('loc_trans: i_prim arg. out of range',at%i_beg(ni)+at%n_hnab(ni)-1,maxatomsproc*at%mx_nab)
          end if
          tpos = at%i_nd_beg(ni)
          nd2 = ahalo%ndimj(ni)
          do nn=1,at%n_hnab(ni)
             i=at%i_prim(at%i_beg(ni)+nn-1)
             istart=ahalo%i_h2d((i-1)*ahalo%ni_in_halo+ni)
             !if(istart>maxatomsproc*at%mx_nab) then
             !   call cq_abort('loc_trans: a arg. out of range',istart,maxatomsproc*at%mx_nab)
             !end if
             nd1 = ahalo%ndimi(i)
             if(istart+nd1*nd2-1>lena) call cq_abort('loc_trans: a arg. out of range',istart+nd1*nd2-1,lena)
             if(tpos+nd1*nd2-1>lenat) call cq_abort('loc_trans: at arg. out of range',tpos+nd1*nd2-1,lenat)
             if(istart<=0) call cq_abort('loc_trans: a arg. out of range',istart)
             if(tpos<=0) call cq_abort('loc_trans: at arg. out of range',tpos)
             do n1=1,nd1
                do n2=1,nd2
                   atrans((n1-1)*nd2+n2-1+tpos)=a((n2-1)*nd1+n1-1+istart)
                end do
             end do
             tpos = tpos+nd2*nd1
          end do
       end do
    elseif(invdir.eq.1) then
       do ni=1,ahalo%ni_in_halo
          if(at%n_hnab(ni)<=0) then
             call cq_abort('loc_trans: no. of neighbours < 1')
          end if
          if(at%i_beg(ni)+at%n_hnab(ni)-1>maxatomsproc*at%mx_nab) then
             call cq_abort('error loc_trans: i_prim arg. out of range',at%i_beg(ni)+at%n_hnab(ni)-1,maxatomsproc*at%mx_nab)
          end if
          tpos = at%i_nd_beg(ni)
          nd2 = ahalo%ndimj(ni)
          do nn=1,at%n_hnab(ni)
             i=at%i_prim(at%i_beg(ni)+nn-1)
             nd1 = ahalo%ndimi(i)
             istart=ahalo%i_h2d((i-1)*ahalo%ni_in_halo+ni)
             !if(istart>maxatomsproc*at%mx_nab) then
             !   call cq_abort('loc_trans: a arg. out of range',istart,maxatomsproc*at%mx_nab)
             !end if
             if(istart+nd1*nd2-1>lena) call cq_abort('loc_trans: a arg. out of range',istart+nd1*nd2-1,lena)
             if(tpos+nd1*nd2-1>lenat) call cq_abort('loc_trans: at arg. out of range',tpos+nd1*nd2-1,lenat)
             if(istart<=0) call cq_abort('loc_trans: a arg. out of range',istart)
             if(tpos<=0) call cq_abort('loc_trans: at arg. out of range',tpos)
             do n1=1,nd1
                do n2=1,nd2
                   a((n2-1)*nd1+n1-1+istart)=atrans((n1-1)*nd2+n2-1+tpos)
                end do
             end do
             tpos = tpos+nd2*nd1
          end do
       end do
    else
       call cq_abort('loc_trans: invdir out of range ',invdir)
    end if
    return
  end subroutine loc_trans
  !!***


  !!****f* multiply_module/prefetch *
  !!
  !!  NAME
  !!   prefetch
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Chooses whether to fetch locally or remotely. If we're hiding
  !!   comms behind calculations, this is called before the matrix mult
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
  !!   2018/10/04 17:31 dave
  !!    Adding tag for MPI compliance
  !!  SOURCE
  !!
  subroutine prefetch(this_part,ahalo,a_b_c,bmat,icall,&
       n_cont,bind_rem,bind,b_rem,lenb_rem,b,myid,ilen2,mx_mpp, &
       parts,prim,gcs,tag)

    ! Module usage
    use datatypes
    use matrix_module
    use basic_types
    use matrix_comms_module
    use comms_module

    implicit none

    ! Passed variables
    ! First, the small group variables
    type(group_set) :: parts
    type(primary_set) :: prim
    type(cover_set) :: gcs
    integer :: mx_mpp
    integer :: this_part,icall,n_cont,myid,ilen2
    type(matrix), dimension(:) :: bmat
    type(matrix_halo) :: ahalo
    type(comms_data) :: a_b_c
    integer(integ), dimension(:)  :: bind_rem,bind
    integer :: lenb_rem, tag
    real(double), dimension(lenb_rem) :: b_rem
    real(double) :: b(:)
    ! Local variables
    integer :: ncover_yz,ind_part,iskip,ind_last
    integer :: inode,ipart,nnode

    ind_part = ahalo%lab_hcell(this_part)
    n_cont=parts%nm_group(ind_part)
    ipart = parts%i_cc2seq(ind_part)
    inode = parts%i_cc2node(ind_part)
    nnode = a_b_c%neigh_node_list(this_part)
    if(inode.eq.myid+1) then ! If this is local, then copy
       icall = 0
       ncover_yz=gcs%ncovery*gcs%ncoverz
       ilen2 = bmat(ipart)%part_nabs
       call Mquest_get_local(ipart,&
            bind_rem,b_rem,lenb_rem,bind,bmat,&
            ind_part,b,myid)
    end if
    if(icall.eq.1) then ! Else fetch the data
       ilen2 = a_b_c%ilen2rec(ipart,nnode)
       call Mquest_get( prim%mx_ngonn, &
            a_b_c%ilen2rec(ipart,nnode),&
            a_b_c%ilen3rec(ipart,nnode),&
            n_cont,inode,ipart,myid,&
            bind_rem,b_rem,lenb_rem,bind,b,&
            a_b_c%istart(ipart,nnode), &
            bmat(1)%mx_abs,parts%mx_mem_grp,tag)
    end if
    return
  end subroutine prefetch
  !!***


  !!****f* multiply_module/check_mkm *
  !!
  !!  NAME
  !!   check_mkm
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Checks dimensions used by m_kern_max/min
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
  !!  SOURCE
  !!
  subroutine check_mkm(kpart,nah_part,kseq,kbeg,ibind_rem,&
     ibpart_rem,nbnab_rem,ibseq_rem,jchbeg,k_off,icall,&
     mx_apart,mx_pcover,mx_ahalo,mx_part,mx_babs,mx_icover)

    use datatypes

    implicit none

    integer :: mx_apart,mx_pcover,mx_icover
    integer :: mx_ahalo,mx_part,mx_babs
    integer :: nbindbeg,kk,jpart,j,jseq,k_off
    integer :: np_in_ahalo,kpart,k,k_in_part,icall
    integer :: kbeg(mx_apart),nah_part(mx_apart)
    integer :: kseq(mx_ahalo)
    integer :: jchbeg(mx_pcover)
    integer(integ) :: ibind_rem(mx_part)
    integer(integ) :: ibpart_rem(mx_part*mx_babs),nbnab_rem(mx_part)
    integer(integ) :: ibseq_rem(mx_part*mx_babs)

    do k=1,nah_part(kpart)
      k_in_part=kseq(kbeg(kpart)+k-1)
      nbindbeg=ibind_rem(k_in_part)
      do j=1,nbnab_rem(k_in_part)
        jpart=ibpart_rem(nbindbeg+j-1)+k_off
        if(jpart.gt.mx_pcover) then
          write(io_lun,911) jpart,mx_pcover
911       format(//'error mat_mult: jchbeg arg. out of range',i12,i6)
          icall = 0
          return
        end if
        jseq=ibseq_rem(nbindbeg+j-1)
        if(jchbeg(jpart)+jseq-1.gt.mx_icover) then
          write(io_lun,912) jchbeg(jpart)+jseq-1,mx_icover
912       format(//'error mat_mult: j_halo arg. out of range',i10,i6)
          icall = 0
          return
        end if
      end do
    end do
    if(icall==0) write(io_lun,*) 'ERROR!!!!!!!!!!!!!!!!!!!!!!!!!'
    return
  end subroutine check_mkm
  !!***


  !!****f* multiply_module/mat_tran *
  !!
  !!  NAME
  !!   mat_tran
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Fetches local transpose from remote nodes
  !!   and puts in place on my node. The nodes from which data are fetched
  !!   are the halo-nodes of my node. I loop over my halo-nodes,
  !!   and from each one I fetch all the local-transpose data in one go.
  !!   I then pass through all these data and use the index
  !!   table itran_addr() to send it to the right address. The index
  !!   itran_addr() is passed to mat_tran in the argument list, and
  !!   should have been created by sbrt ind_tran, which is called
  !!   once only each time atoms are moved.
  !!  INPUTS
  !!   itran_addr(nt): address on my node to which to start putting
  !!    global-transpose data for a given atom pair number nn fetched
  !!    from halo-node nr; in practice, nt=i_pair_indaddr(nr)+nn-1
  !!  USES
  !!
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   25/11/99
  !!  MODIFICATION HISTORY
  !!   20/06/2001 dave
  !!    Added ROBODoc header and cq_abort
  !!   2006/09/07 11:49 dave
  !!    Bug fix for variable function numbers: length now found
  !!    explicitly from pairs%submat
  !!  SOURCE
  !!
  subroutine mat_tran(mx_iprim,mynode,data,atrans,ltrans,pairs,pairsize)

    use datatypes
    use matrix_module
    use trans_module
    use comms_module
    use GenComms, ONLY: cq_abort

    implicit none

    ! Passed variables
    integer(integ) :: mx_iprim
    integer(integ) :: mynode,pairsize
    real(double) :: ltrans(:)
    real(double) :: atrans(:)
    type(trans_remote) :: data
    type(pair_data) :: pairs(:)

    ! Local variables
    real(double), dimension(:), allocatable :: a_store
    integer(integ) :: nr,nnd_rem,tag,n,ic,nn,n1,n2,stat,ierr
    integer :: posn, count,lenat,lensub

    lenat = size(atrans)
    do nr=1,data%n_rem_node ! Loop over halo nodes
       lensub = 0
       do nn=1,data%n_pair(nr)
          lensub = lensub + pairs(nr)%submat(nn)
       end do
       allocate(a_store(lensub),STAT=stat)
       if(stat/=0) call cq_abort("Error allocating a_store in mat_tran ",pairs(nr)%length)
       call fetch_trans_data(a_store,ltrans,data,lensub,nr,mynode,mx_iprim)
       ic=data%i_pair_addr(nr)
       posn = 0
       if(data%n_pair(nr)>data%mx_pair) call cq_abort('Pair error: ',data%n_pair(nr),data%mx_pair)
       do nn=1,data%n_pair(nr)   ! Put it in place
          count = pairs(nr)%submat(nn)
          if(data%itran_addr(ic+nn-1)+count-1>lenat) then
             write(*,*) nr,ic,nn,count,data%itran_addr(ic+nn-1),lenat
             call cq_abort('mat_tran: itran_addr wrong ', data%itran_addr(ic+nn-1)+count-1,lenat)
          end if
          atrans(data%itran_addr(ic+nn-1):data%itran_addr(ic+nn-1)+count-1)=a_store(posn+1:posn+count)
          posn = posn+count
       end do
       deallocate(a_store)
    end do
    return
  end subroutine mat_tran
  !!***


  !!****f* multiply_module/matrix_add *
  !!
  !!  NAME
  !!   matrix_add - adds different range matrices
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Adds two different range matrices and returns the
  !!   result in the first, so A = alpha*A + beta*B
  !!
  !!   Changes the order of arguments passed to addscan so
  !!   that the shorter range matrix is always first - but
  !!   in such a way that the result is put into A (so if
  !!   A is first, dr = 0, and if A is second, dr = 1)
  !!  INPUTS
  !!
  !!
  !!  USES
  !!   addscan
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   25/11/99
  !!  MODIFICATION HISTORY
  !!   20/06/2001 dave
  !!    Added ROBODoc header
  !!  SOURCE
  !!
  subroutine matrix_add(alpha,a,lena,amat,beta,b,lenb,bmat,prim,myid)

    use datatypes
    use matrix_module,ONLY:  matrix
    use basic_types, ONLY: primary_set

    implicit none

    ! Passed variables
    type(matrix) :: amat(:),bmat(:)
    type(primary_set) :: prim
    integer :: lena, lenb
    real(double) :: a(lena),b(lenb)
    real(double) :: alpha,beta
    integer      :: myid

    ! Local variables
    integer :: nn,j,nnd

    nnd = myid+1
    do nn=1,prim%groups_on_node
       if(prim%nm_nodgroup(nn)>0) then
          do j=1,prim%nm_nodgroup(nn)
             if(amat(nn)%n_nab(j)<=bmat(nn)%n_nab(j)) then
                call addscan(prim,nn,j,amat,bmat,a,lena,b,lenb,alpha,beta,myid,0)
             else
                call addscan(prim,nn,j,bmat,amat,b,lenb,a,lena,beta,alpha,myid,1)
             end if
          end do ! j=1,nc_nodpart
       end if ! if(nc_nodpart(nn)>0)
    end do ! nn=1,np_on_node
  end subroutine matrix_add
  !!***


  !!****f* multiply_module/addscan *
  !!
  !!  NAME
  !!   addscan
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Scans through the neighbours of A and B and either adds b to a (dr==0)
  !!   of a to b (dr==1).  Assumes that A is shorter than B.
  !!
  !!   N.B. This only works for the jth atom in the nnth partition of this
  !!   processor's primary set - it's called many times by matrix_add
  !!
  !!   The algorithm is as follows.  Start at the beginning of both matrices'
  !!   lists, and increment B list until it's greater than or equal to A list.
  !!   If equal, check ALL indices, and if they agree, add and continue on.
  !!   If B list is greater than A, then increment A and continue on.
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
  !!  SOURCE
  !!
  subroutine addscan(prim,nn,j,amat,bmat,a,lena,b,lenb,alpha,beta,myid,dr)

    use datatypes
    use matrix_module
    use basic_types
    use GenComms, ONLY: cq_abort

    implicit none

    ! Passed variables
    type(matrix) :: amat(:),bmat(:)
    type(primary_set) :: prim
    integer :: lena, lenb
    real(double) :: a(lena),b(lenb)
    real(double) :: alpha,beta
    integer      :: nn,j,myid,aposn,bposn,dr

    ! Local variables
    integer :: marker,ni,nb,nnd,apos_arr,nd1,nd2, bsize

    nnd = myid+1
    marker = bmat(nn)%i_acc(j)  ! Indexes b
    aposn=amat(nn)%nd_offset+amat(nn)%i_nd_acc(j)  ! Posn in a
    bposn=bmat(nn)%nd_offset+bmat(nn)%i_nd_acc(j)  ! Posn in b
    ni = prim%nm_nodbeg(nn)+j-1  ! Primary atom
    nd1 = amat(nn)%ndimi(j)
    if(amat(nn)%n_nab(j)>0) then  ! Scan over neighbours of j
       do nb=1,amat(nn)%n_nab(j)
          apos_arr = amat(nn)%i_acc(j)+nb-1
          nd2=amat(nn)%ndimj(apos_arr)
          bsize = nd1*nd2
          ! Increment B list until not greater than A
          do while(bmat(nn)%i_part(marker)<amat(nn)%i_part(apos_arr))
             bposn = bposn + nd1*bmat(nn)%ndimj(marker)
             marker = marker+1
          end do
          ! Check against A
          if(bmat(nn)%i_part(marker)>amat(nn)%i_part(apos_arr)) then
             cycle ! We've gone past present neighbour - increment A
          else if(bmat(nn)%i_part(marker)==amat(nn)%i_part(apos_arr)) then
             ! Check through seq and part nos
             do while(bmat(nn)%i_seq(marker)<amat(nn)%i_seq(apos_arr).AND.&
                  bmat(nn)%i_part(marker)==amat(nn)%i_part(apos_arr))
                bposn = bposn + nd1*bmat(nn)%ndimj(marker)
                marker = marker+1
             end do
             ! If either partition no or seq no is too big, increment A
             if(bmat(nn)%i_part(marker)>amat(nn)%i_part(apos_arr)) cycle
             if(bmat(nn)%i_seq(marker)>amat(nn)%i_seq(apos_arr)) then
                cycle
             else if(bmat(nn)%i_part(marker)==amat(nn)%i_part(apos_arr) &
                  .AND.bmat(nn)%i_seq(marker)==amat(nn)%i_seq(apos_arr)) &
                  then
                ! We have found a match - decide which matrix to store result in
                if(dr==0) then
                   a(aposn:aposn+bsize-1)= &
                        alpha*a(aposn:aposn+bsize-1)+beta*b(bposn:bposn+bsize-1)
                else if(dr==1) then
                   b(bposn:bposn+bsize-1)= &
                        beta*b(bposn:bposn+bsize-1)+alpha*a(aposn:aposn+bsize-1)
                end if
                marker=marker+1
                bposn = bposn + bsize
                if(bposn-1>lenb) call cq_abort("Length error for B in addscan: ",bposn-1,lenb)
             end if
          end if
          aposn = aposn+bsize
          if(aposn-1>lena) call cq_abort("Length error for A in addscan: ",aposn-1,lena)
       end do ! nb=1,amat(nn)%n_nab(j)
    end if ! if(amat(nn)%n_nab(j)>0)
  end subroutine addscan
  !!***

end module multiply_module
