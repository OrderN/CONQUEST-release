! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
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
!!  SOURCE
!!
module multiply_module

  use global_module, ONLY: io_lun
  use timer_stdclocks_module, ONLY: start_timer,stop_timer,tmr_std_allocation

  implicit none

  ! RCS tag for object file identification 
  character(len=80), save, private :: RCSid = "$Id$"
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
!!    Changed dimension of nreqs to be mx_nponn*mx_neighbour_nodes (the actual maximum that will be required) and
!!    converted atrans to being dynamically allocated rather than automatic.
!!   10:44, 06/03/2003 drb 
!!    Changed size of nreqs (temporary fix)
!!   09:19, 11/05/2005 dave 
!!    Removed unnecessary MPI_Wait (second one on nreqs)
!!   2008/05/22 ast
!!    Added timers
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
    use GenComms, ONLY: my_barrier, cq_abort, mtime
    use mpi

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
    integer :: icall2,stat,ilen2,lenb_rem
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

    logical flag,call_flag
    real(double) :: t0,t1

    call start_timer(tmr_std_allocation)
    if(iprint_mat>3.AND.myid==0) t0 = mtime()    
    ! Allocate memory for the elements
    allocate(ibpart_rem(a_b_c%parts%mx_mem_grp*a_b_c%bmat(1)%mx_abs),STAT=stat)
    if(stat/=0) call cq_abort('mat_mult: error allocating ibpart_rem')
    !allocate(atrans(a_b_c%amat(1)%length),STAT=stat)
    allocate(atrans(lena),STAT=stat)
    if(stat/=0) call cq_abort('mat_mult: error allocating atrans')
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
    endif
    ! Start the send procedure
    !write(io_lun,*) 'Sending'
    call Mquest_start_send(a_b_c,b,nreqs,myid,a_b_c%prim%mx_ngonn,sends)
    !write(io_lun,*) 'Returned ',a_b_c%ahalo%np_in_halo,myid
    ncover_yz=a_b_c%gcs%ncovery*a_b_c%gcs%ncoverz
    do kpart = 1,a_b_c%ahalo%np_in_halo  ! Main loop
       !write(io_lun,*) 'Part: ',kpart,myid
       icall=1
       ind_part = a_b_c%ahalo%lab_hcell(kpart)
       !write(io_lun,*) 'ind_part: ',ind_part
       if(kpart>1) then  ! Is it a periodic image of the previous partition ?
          if(ind_part.eq.a_b_c%ahalo%lab_hcell(kpart-1)) then 
             icall=0
          else ! Get the data
             !write(io_lun,*) myid,' seq: ',size(a_b_c%parts%i_cc2seq)
             ipart = a_b_c%parts%i_cc2seq(ind_part)
             !write(io_lun,*) myid,' Alloc b_rem part: ',ipart
             nnode = a_b_c%comms%neigh_node_list(kpart)
             !write(io_lun,*) myid,' Alloc b_rem node: ',nnode
             !write(io_lun,*) myid,' Alloc b_rem icc: ',a_b_c%parts%i_cc2node(ind_part)
             !write(io_lun,*) myid,' Alloc b_rem alloc: ',allocated(b_rem)
             if(allocated(b_rem)) deallocate(b_rem)
             if(a_b_c%parts%i_cc2node(ind_part)==myid+1) then
                lenb_rem = a_b_c%bmat(ipart)%part_nd_nabs
             else
                lenb_rem = a_b_c%comms%ilen3rec(ipart,nnode)
             end if
             allocate(b_rem(lenb_rem))
             call prefetch(kpart,a_b_c%ahalo,a_b_c%comms,a_b_c%bmat,icall,&  
                  n_cont,part_array,a_b_c%bindex,b_rem,lenb_rem,b,myid,ilen2,&
                  mx_msg_per_part,a_b_c%parts,a_b_c%prim,a_b_c%gcs)
             !write(io_lun,*) 'b_rem: ',lenb_rem
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
             endif
             ! Create ibpart_rem
             call end_part_comms(myid,n_cont,nbnab_rem,ibind_rem,npxyz_rem,&
                  ibpart_rem,ncover_yz,a_b_c%gcs%ncoverz)
          endif
       else ! Get the data
          !write(io_lun,*) myid,' seq: ',size(a_b_c%parts%i_cc2seq)
          ipart = a_b_c%parts%i_cc2seq(ind_part)
          !write(io_lun,*) myid,' Alloc b_rem part: ',ipart
          nnode = a_b_c%comms%neigh_node_list(kpart)
          !write(io_lun,*) myid,' Alloc b_rem node: ',nnode
          !write(io_lun,*) myid,' Alloc b_rem icc: ',a_b_c%parts%i_cc2node(ind_part)
          !write(io_lun,*) myid,' Alloc b_rem alloc: ',allocated(b_rem)
          if(allocated(b_rem)) deallocate(b_rem)
          if(a_b_c%parts%i_cc2node(ind_part)==myid+1) then
             lenb_rem = a_b_c%bmat(ipart)%part_nd_nabs
          else
             lenb_rem = a_b_c%comms%ilen3rec(ipart,nnode)
          end if
          call start_timer(tmr_std_allocation)
          allocate(b_rem(lenb_rem))
          call stop_timer(tmr_std_allocation)
          call prefetch(kpart,a_b_c%ahalo,a_b_c%comms,a_b_c%bmat,icall,& 
               n_cont,part_array,a_b_c%bindex,b_rem,lenb_rem,b,myid,ilen2,&
               mx_msg_per_part,a_b_c%parts,a_b_c%prim,a_b_c%gcs)
          lenb_rem = size(b_rem)
          !write(io_lun,*) 'b_rem: ',lenb_rem
          ! Now point the _rem variables at the appropriate parts of the array 
          ! where we will receive the data
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
             call cq_abort('Error pointing to part_array !',kpart)
          endif
          call end_part_comms(myid,n_cont,nbnab_rem,ibind_rem,npxyz_rem,&
               ibpart_rem,ncover_yz,a_b_c%gcs%ncoverz)
       endif ! End of the "if this isn't the first partition" loop
       k_off=a_b_c%ahalo%lab_hcover(kpart) ! --- offset for pbcs
       icall2=1
       ! Check dimensions to be used in m_kern_min/max 
       !call check_mkm(kpart,a_b_c%ahalo%nh_part,a_b_c%ahalo%j_seq, &
       !     a_b_c%ahalo%j_beg,&
       !     ibind_rem,ibpart_rem,nbnab_rem,&
       !     ibseq_rem,a_b_c%chalo%i_hbeg,k_off,icall2,&
       !     a_b_c%ahalo%mx_part,a_b_c%gcs%mx_gcover,a_b_c%ahalo%mx_halo, &
       !     a_b_c%parts%mx_mem_grp, &
       !     a_b_c%bmat(1)%mx_abs,a_b_c%gcs%mx_mcover)
       !if(icall2.eq.1) then  ! If check is OK, do the mult
       if(UseGemm) then
          if(a_b_c%mult_type.eq.1) then  ! C is full mult
             call m_kern_maxGEMM( k_off,kpart,ib_nd_acc_rem, ibind_rem,nbnab_rem,&
                  ibpart_rem,ibseq_rem,ibndimj_rem,&
                  atrans,b_rem,c,a_b_c%ahalo,a_b_c%chalo,a_b_c%ltrans,&
                  a_b_c%bmat(1)%mx_abs,a_b_c%parts%mx_mem_grp, &
                  a_b_c%prim%mx_iprim, lena, lenb_rem, lenc)
          else if(a_b_c%mult_type.eq.2) then ! A is partial mult
             call m_kern_minGEMM( k_off,kpart,ib_nd_acc_rem, ibind_rem,nbnab_rem,&
                  ibpart_rem,ibseq_rem,ibndimj_rem,&
                  atrans,b_rem,c,a_b_c%ahalo,a_b_c%chalo,a_b_c%ltrans,&
                  a_b_c%bmat(1)%mx_abs,a_b_c%parts%mx_mem_grp, &
                  a_b_c%prim%mx_iprim, lena, lenb_rem, lenc)
          endif
       else
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
          endif
       end if
       !else
       !   call cq_abort('mat_mult: error in check_mkm ',kpart)
       !endif
    enddo ! End of the kpart=1,ahalo%np_in_halo loop !
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
    endif
    call my_barrier
    call start_timer(tmr_std_allocation)
    deallocate(atrans,STAT=stat)
    if(stat/=0) call cq_abort('mat_mult: error deallocating atrans')
    call stop_timer(tmr_std_allocation)
    call my_barrier
    call start_timer(tmr_std_allocation)
    deallocate(ibpart_rem,STAT=stat)
    if(stat/=0) call cq_abort('mat_mult: error deallocating ibpart_rem')
    call stop_timer(tmr_std_allocation)
    call my_barrier
    !deallocate(b_rem,STAT=stat)
    !if(stat/=0) call cq_abort('mat_mult: error deallocating b_rem')
    !call my_barrier
    if(iprint_mat>3.AND.myid==0) then
       t1 = mtime()
       write(io_lun,*) 'mult time: ',t1-t0
    end if
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
!    enddo
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
          endif
          if(at%i_beg(ni)+at%n_hnab(ni)-1>maxatomsproc*at%mx_nab) then
             call cq_abort('loc_trans: i_prim arg. out of range',at%i_beg(ni)+at%n_hnab(ni)-1,maxatomsproc*at%mx_nab)
          endif
          tpos = at%i_nd_beg(ni)
          nd2 = ahalo%ndimj(ni)
          do nn=1,at%n_hnab(ni)
             i=at%i_prim(at%i_beg(ni)+nn-1)
             istart=ahalo%i_h2d((i-1)*ahalo%ni_in_halo+ni)
             !if(istart>maxatomsproc*at%mx_nab) then
             !   call cq_abort('loc_trans: a arg. out of range',istart,maxatomsproc*at%mx_nab)
             !endif
             nd1 = ahalo%ndimi(i)
             if(istart+nd1*nd2-1>lena) call cq_abort('loc_trans: a arg. out of range',istart+nd1*nd2-1,lena)
             if(tpos+nd1*nd2-1>lenat) call cq_abort('loc_trans: at arg. out of range',tpos+nd1*nd2-1,lenat)
             if(istart<=0) call cq_abort('loc_trans: a arg. out of range',istart)
             if(tpos<=0) call cq_abort('loc_trans: at arg. out of range',tpos)
             do n1=1,nd1
                do n2=1,nd2
                   atrans((n1-1)*nd2+n2-1+tpos)=a((n2-1)*nd1+n1-1+istart)
                enddo
             enddo
             tpos = tpos+nd2*nd1
          enddo
       enddo
    elseif(invdir.eq.1) then
       do ni=1,ahalo%ni_in_halo
          if(at%n_hnab(ni)<=0) then
             call cq_abort('loc_trans: no. of neighbours < 1')
          endif
          if(at%i_beg(ni)+at%n_hnab(ni)-1>maxatomsproc*at%mx_nab) then
             call cq_abort('error loc_trans: i_prim arg. out of range',at%i_beg(ni)+at%n_hnab(ni)-1,maxatomsproc*at%mx_nab)
          endif
          tpos = at%i_nd_beg(ni)
          nd2 = ahalo%ndimj(ni)
          do nn=1,at%n_hnab(ni)
             i=at%i_prim(at%i_beg(ni)+nn-1)
             nd1 = ahalo%ndimi(i)
             istart=ahalo%i_h2d((i-1)*ahalo%ni_in_halo+ni)
             !if(istart>maxatomsproc*at%mx_nab) then
             !   call cq_abort('loc_trans: a arg. out of range',istart,maxatomsproc*at%mx_nab)
             !endif
             if(istart+nd1*nd2-1>lena) call cq_abort('loc_trans: a arg. out of range',istart+nd1*nd2-1,lena)
             if(tpos+nd1*nd2-1>lenat) call cq_abort('loc_trans: at arg. out of range',tpos+nd1*nd2-1,lenat)
             if(istart<=0) call cq_abort('loc_trans: a arg. out of range',istart)
             if(tpos<=0) call cq_abort('loc_trans: at arg. out of range',tpos)
             do n1=1,nd1
                do n2=1,nd2
                   a((n2-1)*nd1+n1-1+istart)=atrans((n1-1)*nd2+n2-1+tpos)
                enddo
             enddo
             tpos = tpos+nd2*nd1
          enddo
       enddo
    else
       call cq_abort('loc_trans: invdir out of range ',invdir)
    endif
    return
  end subroutine loc_trans
!!***

!!****f* multiply_module/m_kern_maxGEMM *
!!
!!  NAME 
!!   m_kern_max
!!  USAGE
!! 
!!  PURPOSE
!!   multiplication kernel for maximal case and reductions thereof.
!!
!!   nah_part:       the number of atoms in the A-halo that are contained
!!                   in the current partition K.
!!
!!   kseq(k):        the unpruned sequence number of atom k in the current
!!                   partition. `Unpruned' means that the sequence
!!                   number counts all atoms in the partition, and not just
!!                   those in the A-halo. The atom in question is specified
!!                   by its `pruned' sequence number k.
!! 
!!   k_in_part:      temporary variable for current value of kseq(k).
!! 
!!   k_halo(.):      the A-halo sequence number of atom k. The latter is
!!                   specified by its partition number and its pruned
!!                   sequence number k in that partition.
!! 
!!   k_in_halo:      temporary variable for current value of k_halo(.)
!! 
!!   kbeg(kpart):    specifies address in k_halo(.) for start of information
!!                   about atoms in partition kpart. The label kpart
!!                   goes over all partitions containing atoms in the
!!                   A-halo.
!!
!!   kpart:          A-halo seq. no. of current partition K
!!
!!   nahnab(k_in_halo):
!!                   number of atoms in primary set that are A-neighbours of
!!                   a given atom in the A-halo, the latter being given by its
!!                   A-halo sequence number k_in_halo.
!!
!!   i_prim(.):      sequence number of atom i in primary set, this
!!                   atom being specified as ith neighbour of atom
!!                   k in the A-halo
!!
!!   ibeg(k_in_halo):
!!                   address in i_prim(.) where index data for atom k
!!                   start, the latter being specified by k_in_halo.
!!
!!   iaaddr(k_in_halo):
!!                   address in array a where data for atom k start,
!!                   the latter being specified by k_in_halo.
!!
!!   nbnab(k_in_part):
!!                   number of B-neighbours of atom k, the latter being
!!                   specified by its (unpruned) seq. no. in partition K.
!!
!!   ibaddr(k_in_part):
!!                   address in array b where data for atom k start.
!!
!!   jch2cad(.):     address in array c where data for (i,j) pair start.
!!
!!   icad:           address in jch2cad where info for atom i starts.
!!
!!   jbnab2ch(.):    the C-halo seq. no. of atom j, the latter being
!!                   specified as the jth B-neighbour of atom k.
!!
!!   ni_in_chalo:        total number of atoms in the C-halo.
!!
!!   ibpart(.):      partition to which each atom j belongs, the latter
!!                   being specified as jth B-neighbour of atom k.
!!
!!   ibseq(.):       unpruned sequence number in partition for each
!!                   atom j, the latter being specified as in ibpart.
!!
!!   ibindaddr(k_in_part):
!!                   address in arrays ibpart(.) and ibseq(.) where
!!                   data for atom k start.
!!
!!   k_off:          offset in partition labelling to account for p.b.c.
!!
!!   j_halo(.):      C-halo sequence number of atom j, the latter being
!!                   specified in partition labelling.
!!
!!   jchbeg(jpart):   address in array j_halo(.) where data for
!!                   partition jpart start.
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   M.J.Gillan/D.R.Bowler
!!  CREATION DATE
!!   25/11/99
!!  MODIFICATION HISTORY
!!   20/06/2001 dave
!!    Added ROBODoc header
!!  SOURCE
!!
  subroutine m_kern_maxGEMM(k_off,kpart,&
       ib_nd_acc,ibaddr,nbnab,ibpart,ibseq,bndim2,a,b,c,ahalo,chalo,at,&
       mx_absb,mx_part,mx_iprim, lena, lenb, lenc,debug)

    use datatypes
    use matrix_module
    use basic_types, ONLY: primary_set
    use primary_module, ONLY: bundle
    use numbers, ONLY: zero

    implicit none

    ! Passed variables
    type(matrix_halo)::ahalo,chalo
    type(matrix_trans)::at
    integer :: mx_absb,mx_part,mx_iprim,lena, lenb, lenc
    real(double) :: a(lena)
    real(double) :: b(lenb)
    real(double) :: c(lenc)

    integer :: kpart,k_off
    integer, OPTIONAL :: debug
    ! Remote indices
    integer(integ) :: ib_nd_acc(mx_part)
    integer(integ) :: ibaddr(mx_part)
    integer(integ) :: nbnab(mx_part)
    integer(integ) :: ibpart(mx_part*mx_absb)
    integer(integ) :: ibseq(mx_part*mx_absb)
    integer(integ) :: bndim2(mx_part*mx_absb)
    ! Local variables
    integer :: jbnab2ch(mx_absb)  ! Automatic array
    integer :: nbkbeg,k,k_in_part,k_in_halo,j,jpart,jseq
    integer :: i,nabeg,i_in_prim,icad,nbbeg,j_in_halo,ncbeg
    integer :: n1,n2,n3,nb_nd_kbeg
    integer :: nd1,nd2,nd3
    integer :: naaddr,nbaddr,ncaddr
    real(double), allocatable, dimension(:,:) :: tempb,tempa, tempc
    integer :: sofar, maxlen, max2,prend1
    external :: dgemm

    allocate(tempa(1,1),tempc(1,1))
    do k=1,ahalo%nh_part(kpart)  ! Loop over atoms k in current A-halo partn
       k_in_halo=ahalo%j_beg(kpart)+k-1
       k_in_part=ahalo%j_seq(k_in_halo)
       nbkbeg=ibaddr(k_in_part)
       nb_nd_kbeg=ib_nd_acc(k_in_part)
       nd3 = ahalo%ndimj(k_in_halo)
       !if(PRESENT(debug)) write(21+debug,*) 'Details1: ',k,nb_nd_kbeg
       ! transcription of j from partition to C-halo labelling
       maxlen = 0
       do j=1,nbnab(k_in_part)
          jpart=ibpart(nbkbeg+j-1)+k_off
          jseq=ibseq(nbkbeg+j-1)
          jbnab2ch(j)=chalo%i_halo(chalo%i_hbeg(jpart)+jseq-1)
          maxlen = maxlen + bndim2(nbkbeg+j-1)
       enddo
       nabeg=at%i_nd_beg(k_in_halo)  ! nabeg=at%i_beg(k_in_halo)
       allocate(tempb(nd3,maxlen))
       prend1=0
       do i=1,at%n_hnab(k_in_halo)  ! Loop over primary-set A-neighbours of k
          !nabeg=at%i_beg(k_in_halo)+i-1
          i_in_prim=at%i_prim(at%i_beg(k_in_halo)+i-1)
          nd1 = ahalo%ndimi(i_in_prim)
          if(nd1/=prend1) then
             deallocate(tempc,tempa)          
             allocate(tempa(nd1,nd3),tempc(nd1,maxlen))
             !allocate(tempa(nd3,nd1),tempc(nd1,maxlen))
          end if
          tempa = zero
          tempb = zero
          tempc = zero
          do n1=1,nd1
             naaddr=nabeg+nd3*(n1-1)
             do n3 = 1,nd3
                tempa(n1,n3) = a(naaddr+n3-1)
                !tempa(n3,n1) = a(naaddr+n3-1)
             end do
          end do
          icad=(i_in_prim-1)*chalo%ni_in_halo
          nbbeg=nb_nd_kbeg
          sofar = 0
          do j=1,nbnab(k_in_part) ! Loop over B-neighbours of atom k
             !nbbeg=nbkbeg+j-1
             nd2 = bndim2(nbkbeg+j-1)
             j_in_halo=jbnab2ch(j)
             if(j_in_halo.ne.0) then
                ncbeg=chalo%i_h2d(icad+j_in_halo)
                !nd2 = chalo%ndimj(j_in_halo)
                if(ncbeg/=0) then  ! multiplication of ndim x ndim blocks 
                   !if (PRESENT(debug)) write(21+debug,*) &
                   !     'Details2: ',j,nd2,(nabeg-1)/(nd1*nd3),(ncbeg-1)/(nd1*nd2),(nbbeg-1)/(nd2*nd3)
!DIR$ NOPATTERN
                   !!  do n2=1,nd2
                   !!     nbaddr = nbbeg+nd3*(n2-1)
                   !!     ncaddr = ncbeg+nd1*(n2-1)
                   !!     do n1=1,nd1
                   !!        naaddr=nabeg+nd3*(n1-1)
                   !!        do n3=1,nd3
                   !!           c(ncaddr+n1-1) = c(ncaddr+n1-1) &
                   !!              +a(naaddr+n3-1)*b(nbaddr+n3-1)
                   !!        enddo
                   !!     enddo
                   !!  enddo
                   do n2=1,nd2
                      nbaddr = nbbeg+nd3*(n2-1)
                      do n3 = 1,nd3
                         tempb(n3,sofar+n2) = b(nbaddr+n3-1)
                      enddo
                   end do
                   sofar = sofar + nd2
                endif
             endif ! End of if(j_in_halo.ne.0)
             nbbeg = nbbeg + nd3*nd2
          enddo ! End of 1,nbnab
          if(sofar>0) then
             ! m,n,k,alpha,a,lda,b,ldb,beta,c,ldc
             !call dgemm('t','n',nd1,sofar,nd3,1.0_double,tempa,nd3,tempb,nd3,0.0_double,tempc,nd1)
             call dgemm('n','n',nd1,sofar,nd3,1.0_double,tempa,nd1,tempb,nd3,0.0_double,tempc,nd1)
          end if
          sofar = 0
          do j=1,nbnab(k_in_part) ! Loop over B-neighbours of atom k
             nd2 = bndim2(nbkbeg+j-1)
             j_in_halo=jbnab2ch(j)
             if(j_in_halo.ne.0) then
                ncbeg=chalo%i_h2d(icad+j_in_halo)
                if(ncbeg/=0) then  ! multiplication of ndim x ndim blocks 
                   do n2=1,nd2
                      ncaddr = ncbeg+nd1*(n2-1)
                      do n1=1,nd1
                         c(ncaddr+n1-1) = c(ncaddr+n1-1) + tempc(n1,sofar+n2)
                      end do
                   end do
                   sofar = sofar + nd2
                end if
             end if
          end do
          nabeg = nabeg + nd1*nd3
       enddo ! End of 1,at%n_hnab
       deallocate(tempb)          
    enddo ! End of k=1,nahpart
    if(allocated(tempa)) deallocate(tempa)          
    if(allocated(tempc)) deallocate(tempc)          
    return
  end subroutine m_kern_maxGEMM
!!***

!!****f* multiply_module/m_kern_minGEMM *
!!
!!  NAME 
!!   m_kern_min
!!  USAGE
!! 
!!  PURPOSE
!!   multiplication kernel for minimal case and extensions thereof.
!!
!!   nah_part:       the number of atoms in the A-halo that are contained
!!                   in the current partition K.
!!
!!   kseq(k):        the unpruned sequence number of atom k in the current
!!                   partition. `Unpruned' means that the sequence
!!                   number counts all atoms in the partition, and not just
!!                   those in the A-halo. The atom in question is specified
!!                   by its `pruned' sequence number k.
!! 
!!   k_in_part:      temporary variable for current value of kseq(k).
!! 
!!   k_halo(.):      the A-halo sequence number of atom k. The latter is
!!                   specified by its partition number and its pruned
!!                   sequence number k in that partition.
!! 
!!   k_in_halo:      temporary variable for current value of k_halo(.)
!! 
!!   kbeg(kpart):    specifies address in k_halo(.) for start of information
!!                   about atoms in partition kpart. The label kpart
!!                   goes over all partitions containing atoms in the
!!                   A-halo.
!!
!!   kpart:          A-halo seq. no. of current partition K
!!
!!   nahnab(k_in_halo):
!!                   number of atoms in primary set that are A-neighbours of
!!                   a given atom in the A-halo, the latter being given by its
!!                   A-halo sequence number k_in_halo.
!!
!!   i_prim(.):      sequence number of atom i in primary set, this
!!                   atom being specified as ith neighbour of atom
!!                   k in the A-halo
!!
!!   ibeg(k_in_halo):
!!                   address in i_prim(.) where index data for atom k
!!                   start, the latter being specified by k_in_halo.
!!
!!   iaaddr(k_in_halo):
!!                   address in array a where data for atom k start,
!!                   the latter being specified by k_in_halo.
!!
!!   nbnab(k_in_part):
!!                   number of B-neighbours of atom k, the latter being
!!                   specified by its (unpruned) seq. no. in partition K.
!!
!!   ibaddr(k_in_part):
!!                   address in array b where data for atom k start.
!!
!!   jch2cad(.):     address in array ! where data for (i,j) pair start.
!!
!!   icad:           address in jch2cad where info for atom i starts.
!!
!!   jbnab2ch(.):    the C-halo seq. no. of atom j, the latter being
!!                   specified as the jth B-neighbour of atom k.
!!
!!   ni_in_chalo:        total number of atoms in the C-halo.
!!
!!   ibpart(.):      partition to which each atom j belongs, the latter
!!                   being specified as jth B-neighbour of atom k.
!!
!!   ibseq(.):       unpruned sequence number in partition for each
!!                   atom j, the latter being specified as in ibpart.
!!
!!   ibindaddr(k_in_part):
!!                   address in arrays ibpart(.) and ibseq(.) where
!!                   data for atom k start.
!!
!!   k_off:          offset in partition labelling to account for p.b.c.
!!
!!   j_halo(.):      C-halo sequence number of atom j, the latter being
!!                   specified in partition labelling.
!!
!!   jchbeg(jpart):   address in array j_halo(.) where data for
!!                   partition jpart start.
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   M.J.Gillan/D.R.Bowler
!!  CREATION DATE
!!   25/11/99
!!  MODIFICATION HISTORY
!!   20/06/2001 dave
!!    Added ROBODoc header
!!  SOURCE
!!
  subroutine m_kern_minGEMM(k_off,kpart,&
       ib_nd_acc,ibaddr,nbnab,ibpart,ibseq,bndim2,a,b,c,ahalo,chalo,at,&
       mx_absb,mx_part,mx_iprim, lena, lenb, lenc)

    use datatypes
    use matrix_module
    use basic_types, ONLY: primary_set
    use primary_module, ONLY: bundle

    implicit none

    ! Passed variables
    type(matrix_halo)::ahalo,chalo
    type(matrix_trans)::at
    integer :: mx_absb,mx_part,mx_iprim,lena, lenb, lenc
    ! Remember that a is a local transpose
    real(double) :: a(lena)
    real(double) :: b(lenb)
    real(double) :: c(lenc)

    integer :: kpart,k_off
    ! dimension declarations
    integer :: ibaddr(mx_part)
    integer :: ib_nd_acc(mx_part)
    integer :: nbnab(mx_part)
    integer :: ibpart(mx_part*mx_absb)
    integer :: ibseq(mx_part*mx_absb)
    integer :: bndim2(mx_part*mx_absb)
    ! Local variables
    integer :: jbnab2ch(mx_absb)
    integer :: k,k_in_part,k_in_halo,nbkbeg,j,jpart,jseq
    integer :: i,nabeg,i_in_prim,icad,nbbeg,j_in_halo,ncbeg
    integer :: n1,n2,n3,nb_nd_kbeg
    integer :: nd1,nd2,nd3
    integer :: naaddr,nbaddr,ncaddr

    real(double), allocatable, dimension(:,:) :: store!(ndim3,ndim1)
    real(double), allocatable, dimension(:,:) :: tempb,tempc
    integer :: sofar, maxlen
    external :: dgemm

    do k=1,ahalo%nh_part(kpart) ! Loop over atoms k in current A-halo partn
       k_in_halo=ahalo%j_beg(kpart)+k-1
       k_in_part=ahalo%j_seq(k_in_halo)
       nbkbeg=ibaddr(k_in_part)
       nb_nd_kbeg=ib_nd_acc(k_in_part)
       nd3 = ahalo%ndimj(k_in_halo)
       ! transcription of j from partition to C-halo labelling
       maxlen = 0
       do j=1,nbnab(k_in_part)
          jpart=ibpart(nbkbeg+j-1)+k_off
          jseq=ibseq(nbkbeg+j-1)
          jbnab2ch(j)=chalo%i_halo(chalo%i_hbeg(jpart)+jseq-1)
          maxlen = maxlen + bndim2(nbkbeg+j-1)
       enddo
       nabeg=at%i_nd_beg(k_in_halo)  ! nabeg=at%i_beg(k_in_halo)  
       do i=1,at%n_hnab(k_in_halo)  ! Loop over primary-set A-neighbours of k
          !nabeg=at%i_beg(k_in_halo)+i-1
          i_in_prim=at%i_prim(at%i_beg(k_in_halo)+i-1)
          nd1 = ahalo%ndimi(i_in_prim)
          icad=(i_in_prim-1)*chalo%ni_in_halo
          nbbeg=nb_nd_kbeg
          sofar = 0
          allocate(tempb(nd3,maxlen),tempc(maxlen,nd1),store(nd3,nd1))
          do j=1,nbnab(k_in_part)  ! Loop over B-neighbours of atom k
             !nbbeg=nbkbeg+j-1
             nd2 = bndim2(nbkbeg+j-1)
             j_in_halo=jbnab2ch(j)
             if(j_in_halo.ne.0) then
                !nd2 = chalo%ndimj(j_in_halo)
                ncbeg=chalo%i_h2d(icad+j_in_halo)
                if(ncbeg.ne.0) then  ! multiplication of ndim x ndim blocks 
!DIR$ NOPATTERN
                   !do n2=1,nd2
                   !   nbaddr = nbbeg+nd3*(n2-1)
                   !   ncaddr = ncbeg+nd1*(n2-1)
                   !   do n1=1,nd1
                   !      naaddr=nabeg+nd3*(n1-1)
                   !      do n3=1,nd3
                   !         a(naaddr+n3-1) = a(naaddr+n3-1) &
                   !              +c(ncaddr+n1-1)*b(nbaddr+n3-1)
                   !      enddo
                   !   enddo
                   !enddo
                   !!   Alternative
                   do n2=1,nd2
                      nbaddr = nbbeg+nd3*(n2-1)
                      ncaddr = ncbeg+nd1*(n2-1)
                      do n3=1,nd3
                         tempb(n3,sofar+n2) = b(nbaddr+n3-1)
                      end do
                      do n1=1,nd1
                         tempc(sofar+n2,n1) = c(ncaddr+n1-1)
                      end do
                   end do
                   sofar = sofar+nd2
                endif
             endif
             nbbeg = nbbeg + nd3*nd2
          enddo
          if(sofar>0) then
             call dgemm('n','n',nd3,nd1,sofar,1.0_double,tempb,nd3,tempc,maxlen,1.0_double,a(nabeg:),nd3)
          end if
          deallocate(tempb,tempc,store)          
          nabeg = nabeg + nd1*nd3
       enddo
    enddo
    return
  end subroutine m_kern_minGEMM
!!***

!!****f* multiply_module/m_kern_max *
!!
!!  NAME 
!!   m_kern_max
!!  USAGE
!! 
!!  PURPOSE
!!   multiplication kernel for maximal case and reductions thereof.
!!
!!   nah_part:       the number of atoms in the A-halo that are contained
!!                   in the current partition K.
!!
!!   kseq(k):        the unpruned sequence number of atom k in the current
!!                   partition. `Unpruned' means that the sequence
!!                   number counts all atoms in the partition, and not just
!!                   those in the A-halo. The atom in question is specified
!!                   by its `pruned' sequence number k.
!! 
!!   k_in_part:      temporary variable for current value of kseq(k).
!! 
!!   k_halo(.):      the A-halo sequence number of atom k. The latter is
!!                   specified by its partition number and its pruned
!!                   sequence number k in that partition.
!! 
!!   k_in_halo:      temporary variable for current value of k_halo(.)
!! 
!!   kbeg(kpart):    specifies address in k_halo(.) for start of information
!!                   about atoms in partition kpart. The label kpart
!!                   goes over all partitions containing atoms in the
!!                   A-halo.
!!
!!   kpart:          A-halo seq. no. of current partition K
!!
!!   nahnab(k_in_halo):
!!                   number of atoms in primary set that are A-neighbours of
!!                   a given atom in the A-halo, the latter being given by its
!!                   A-halo sequence number k_in_halo.
!!
!!   i_prim(.):      sequence number of atom i in primary set, this
!!                   atom being specified as ith neighbour of atom
!!                   k in the A-halo
!!
!!   ibeg(k_in_halo):
!!                   address in i_prim(.) where index data for atom k
!!                   start, the latter being specified by k_in_halo.
!!
!!   iaaddr(k_in_halo):
!!                   address in array a where data for atom k start,
!!                   the latter being specified by k_in_halo.
!!
!!   nbnab(k_in_part):
!!                   number of B-neighbours of atom k, the latter being
!!                   specified by its (unpruned) seq. no. in partition K.
!!
!!   ibaddr(k_in_part):
!!                   address in array b where data for atom k start.
!!
!!   jch2cad(.):     address in array c where data for (i,j) pair start.
!!
!!   icad:           address in jch2cad where info for atom i starts.
!!
!!   jbnab2ch(.):    the C-halo seq. no. of atom j, the latter being
!!                   specified as the jth B-neighbour of atom k.
!!
!!   ni_in_chalo:        total number of atoms in the C-halo.
!!
!!   ibpart(.):      partition to which each atom j belongs, the latter
!!                   being specified as jth B-neighbour of atom k.
!!
!!   ibseq(.):       unpruned sequence number in partition for each
!!                   atom j, the latter being specified as in ibpart.
!!
!!   ibindaddr(k_in_part):
!!                   address in arrays ibpart(.) and ibseq(.) where
!!                   data for atom k start.
!!
!!   k_off:          offset in partition labelling to account for p.b.c.
!!
!!   j_halo(.):      C-halo sequence number of atom j, the latter being
!!                   specified in partition labelling.
!!
!!   jchbeg(jpart):   address in array j_halo(.) where data for
!!                   partition jpart start.
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   M.J.Gillan/D.R.Bowler
!!  CREATION DATE
!!   25/11/99
!!  MODIFICATION HISTORY
!!   20/06/2001 dave
!!    Added ROBODoc header
!!  SOURCE
!!
  subroutine m_kern_max(k_off,kpart,&
       ib_nd_acc,ibaddr,nbnab,ibpart,ibseq,bndim2,a,b,c,ahalo,chalo,at,&
       mx_absb,mx_part,mx_iprim, lena, lenb, lenc,debug)

    use datatypes
    use matrix_module
    use basic_types, ONLY: primary_set
    use primary_module, ONLY: bundle

    implicit none

    ! Passed variables
    type(matrix_halo)::ahalo,chalo
    type(matrix_trans)::at
    integer :: mx_absb,mx_part,mx_iprim,lena, lenb, lenc
    real(double) :: a(lena)
    real(double) :: b(lenb)
    real(double) :: c(lenc)

    integer :: kpart,k_off
    integer, OPTIONAL :: debug
    ! Remote indices
    integer(integ) :: ib_nd_acc(mx_part)
    integer(integ) :: ibaddr(mx_part)
    integer(integ) :: nbnab(mx_part)
    integer(integ) :: ibpart(mx_part*mx_absb)
    integer(integ) :: ibseq(mx_part*mx_absb)
    integer(integ) :: bndim2(mx_part*mx_absb)
    ! Local variables
    integer :: jbnab2ch(mx_absb)  ! Automatic array
    integer :: nbkbeg,k,k_in_part,k_in_halo,j,jpart,jseq
    integer :: i,nabeg,i_in_prim,icad,nbbeg,j_in_halo,ncbeg
    integer :: n1,n2,n3,nb_nd_kbeg
    integer :: nd1,nd2,nd3
    integer :: naaddr,nbaddr,ncaddr

    do k=1,ahalo%nh_part(kpart)  ! Loop over atoms k in current A-halo partn
       k_in_halo=ahalo%j_beg(kpart)+k-1
       k_in_part=ahalo%j_seq(k_in_halo)
       nbkbeg=ibaddr(k_in_part)
       nb_nd_kbeg=ib_nd_acc(k_in_part)
       nd3 = ahalo%ndimj(k_in_halo)
       !if(PRESENT(debug)) write(21+debug,*) 'Details1: ',k,nb_nd_kbeg
       ! transcription of j from partition to C-halo labelling
       do j=1,nbnab(k_in_part)
          jpart=ibpart(nbkbeg+j-1)+k_off
          jseq=ibseq(nbkbeg+j-1)
          jbnab2ch(j)=chalo%i_halo(chalo%i_hbeg(jpart)+jseq-1)
       enddo
       nabeg=at%i_nd_beg(k_in_halo)  ! nabeg=at%i_beg(k_in_halo)
       do i=1,at%n_hnab(k_in_halo)  ! Loop over primary-set A-neighbours of k
          !nabeg=at%i_beg(k_in_halo)+i-1
          i_in_prim=at%i_prim(at%i_beg(k_in_halo)+i-1)
          nd1 = ahalo%ndimi(i_in_prim)
          icad=(i_in_prim-1)*chalo%ni_in_halo
          nbbeg=nb_nd_kbeg
          do j=1,nbnab(k_in_part) ! Loop over B-neighbours of atom k
             !nbbeg=nbkbeg+j-1
             nd2 = bndim2(nbkbeg+j-1)
             j_in_halo=jbnab2ch(j)
             if(j_in_halo.ne.0) then
                ncbeg=chalo%i_h2d(icad+j_in_halo)
                !nd2 = chalo%ndimj(j_in_halo)
                if(ncbeg/=0) then  ! multiplication of ndim x ndim blocks 
                   !if (PRESENT(debug)) write(21+debug,*) &
                   !     'Details2: ',j,nd2,(nabeg-1)/(nd1*nd3),(ncbeg-1)/(nd1*nd2),(nbbeg-1)/(nd2*nd3)
!DIR$ NOPATTERN
                   do n2=1,nd2
                      nbaddr = nbbeg+nd3*(n2-1)
                      ncaddr = ncbeg+nd1*(n2-1)
                      do n1=1,nd1
                         naaddr=nabeg+nd3*(n1-1)
                         do n3=1,nd3
                            c(ncaddr+n1-1) = c(ncaddr+n1-1) &
                                 +a(naaddr+n3-1)*b(nbaddr+n3-1)
                         enddo
                      enddo
                   enddo
                endif
             endif ! End of if(j_in_halo.ne.0)
             nbbeg = nbbeg + nd3*nd2
          enddo ! End of 1,nbnab
          nabeg = nabeg + nd1*nd3
       enddo ! End of 1,at%n_hnab
    enddo ! End of k=1,nahpart
    return
  end subroutine m_kern_max
!!***

!!****f* multiply_module/m_kern_min *
!!
!!  NAME 
!!   m_kern_min
!!  USAGE
!! 
!!  PURPOSE
!!   multiplication kernel for minimal case and extensions thereof.
!!
!!   nah_part:       the number of atoms in the A-halo that are contained
!!                   in the current partition K.
!!
!!   kseq(k):        the unpruned sequence number of atom k in the current
!!                   partition. `Unpruned' means that the sequence
!!                   number counts all atoms in the partition, and not just
!!                   those in the A-halo. The atom in question is specified
!!                   by its `pruned' sequence number k.
!! 
!!   k_in_part:      temporary variable for current value of kseq(k).
!! 
!!   k_halo(.):      the A-halo sequence number of atom k. The latter is
!!                   specified by its partition number and its pruned
!!                   sequence number k in that partition.
!! 
!!   k_in_halo:      temporary variable for current value of k_halo(.)
!! 
!!   kbeg(kpart):    specifies address in k_halo(.) for start of information
!!                   about atoms in partition kpart. The label kpart
!!                   goes over all partitions containing atoms in the
!!                   A-halo.
!!
!!   kpart:          A-halo seq. no. of current partition K
!!
!!   nahnab(k_in_halo):
!!                   number of atoms in primary set that are A-neighbours of
!!                   a given atom in the A-halo, the latter being given by its
!!                   A-halo sequence number k_in_halo.
!!
!!   i_prim(.):      sequence number of atom i in primary set, this
!!                   atom being specified as ith neighbour of atom
!!                   k in the A-halo
!!
!!   ibeg(k_in_halo):
!!                   address in i_prim(.) where index data for atom k
!!                   start, the latter being specified by k_in_halo.
!!
!!   iaaddr(k_in_halo):
!!                   address in array a where data for atom k start,
!!                   the latter being specified by k_in_halo.
!!
!!   nbnab(k_in_part):
!!                   number of B-neighbours of atom k, the latter being
!!                   specified by its (unpruned) seq. no. in partition K.
!!
!!   ibaddr(k_in_part):
!!                   address in array b where data for atom k start.
!!
!!   jch2cad(.):     address in array ! where data for (i,j) pair start.
!!
!!   icad:           address in jch2cad where info for atom i starts.
!!
!!   jbnab2ch(.):    the C-halo seq. no. of atom j, the latter being
!!                   specified as the jth B-neighbour of atom k.
!!
!!   ni_in_chalo:        total number of atoms in the C-halo.
!!
!!   ibpart(.):      partition to which each atom j belongs, the latter
!!                   being specified as jth B-neighbour of atom k.
!!
!!   ibseq(.):       unpruned sequence number in partition for each
!!                   atom j, the latter being specified as in ibpart.
!!
!!   ibindaddr(k_in_part):
!!                   address in arrays ibpart(.) and ibseq(.) where
!!                   data for atom k start.
!!
!!   k_off:          offset in partition labelling to account for p.b.c.
!!
!!   j_halo(.):      C-halo sequence number of atom j, the latter being
!!                   specified in partition labelling.
!!
!!   jchbeg(jpart):   address in array j_halo(.) where data for
!!                   partition jpart start.
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   M.J.Gillan/D.R.Bowler
!!  CREATION DATE
!!   25/11/99
!!  MODIFICATION HISTORY
!!   20/06/2001 dave
!!    Added ROBODoc header
!!  SOURCE
!!
  subroutine m_kern_min(k_off,kpart,&
       ib_nd_acc,ibaddr,nbnab,ibpart,ibseq,bndim2,a,b,c,ahalo,chalo,at,&
       mx_absb,mx_part,mx_iprim, lena, lenb, lenc)

    use datatypes
    use matrix_module
    use basic_types, ONLY: primary_set
    use primary_module, ONLY: bundle

    implicit none

    ! Passed variables
    type(matrix_halo)::ahalo,chalo
    type(matrix_trans)::at
    integer :: mx_absb,mx_part,mx_iprim,lena, lenb, lenc
    ! Remember that a is a local transpose
    real(double) :: a(lena)
    real(double) :: b(lenb)
    real(double) :: c(lenc)

    integer :: kpart,k_off
    ! dimension declarations
    integer :: ibaddr(mx_part)
    integer :: ib_nd_acc(mx_part)
    integer :: nbnab(mx_part)
    integer :: ibpart(mx_part*mx_absb)
    integer :: ibseq(mx_part*mx_absb)
    integer :: bndim2(mx_part*mx_absb)
    ! Local variables
    integer :: jbnab2ch(mx_absb)
    integer :: k,k_in_part,k_in_halo,nbkbeg,j,jpart,jseq
    integer :: i,nabeg,i_in_prim,icad,nbbeg,j_in_halo,ncbeg
    integer :: n1,n2,n3,nb_nd_kbeg
    integer :: nd1,nd2,nd3
    integer :: naaddr,nbaddr,ncaddr

    do k=1,ahalo%nh_part(kpart) ! Loop over atoms k in current A-halo partn
       k_in_halo=ahalo%j_beg(kpart)+k-1
       k_in_part=ahalo%j_seq(k_in_halo)
       nbkbeg=ibaddr(k_in_part)
       nb_nd_kbeg=ib_nd_acc(k_in_part)
       nd3 = ahalo%ndimj(k_in_halo)
       ! transcription of j from partition to C-halo labelling
       do j=1,nbnab(k_in_part)
          jpart=ibpart(nbkbeg+j-1)+k_off
          jseq=ibseq(nbkbeg+j-1)
          jbnab2ch(j)=chalo%i_halo(chalo%i_hbeg(jpart)+jseq-1)
       enddo
       nabeg=at%i_nd_beg(k_in_halo)  ! nabeg=at%i_beg(k_in_halo)  
       do i=1,at%n_hnab(k_in_halo)  ! Loop over primary-set A-neighbours of k
          !nabeg=at%i_beg(k_in_halo)+i-1
          i_in_prim=at%i_prim(at%i_beg(k_in_halo)+i-1)
          nd1 = ahalo%ndimi(i_in_prim)
          icad=(i_in_prim-1)*chalo%ni_in_halo
          nbbeg=nb_nd_kbeg
          do j=1,nbnab(k_in_part)  ! Loop over B-neighbours of atom k
             !nbbeg=nbkbeg+j-1
             nd2 = bndim2(nbkbeg+j-1)
             j_in_halo=jbnab2ch(j)
             if(j_in_halo.ne.0) then
                !nd2 = chalo%ndimj(j_in_halo)
                ncbeg=chalo%i_h2d(icad+j_in_halo)
                if(ncbeg.ne.0) then  ! multiplication of ndim x ndim blocks 
!DIR$ NOPATTERN
                   do n2=1,nd2
                      nbaddr = nbbeg+nd3*(n2-1)
                      ncaddr = ncbeg+nd1*(n2-1)
                      do n1=1,nd1
                         naaddr=nabeg+nd3*(n1-1)
                         do n3=1,nd3
                            a(naaddr+n3-1) = a(naaddr+n3-1) &
                                 +c(ncaddr+n1-1)*b(nbaddr+n3-1)
                         enddo
                      enddo
                   enddo
                endif
             endif
             nbbeg = nbbeg + nd3*nd2
          enddo
          nabeg = nabeg + nd1*nd3
       enddo
    enddo
    return
  end subroutine m_kern_min
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
!!  SOURCE
!!
  subroutine prefetch(this_part,ahalo,a_b_c,bmat,icall,&
       n_cont,bind_rem,bind,b_rem,lenb_rem,b,myid,ilen2,mx_mpp, &
       parts,prim,gcs)

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
    integer :: lenb_rem
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
    endif
    if(icall.eq.1) then ! Else fetch the data
       ilen2 = a_b_c%ilen2rec(ipart,nnode)
       call Mquest_get( prim%mx_ngonn, &
            a_b_c%ilen2rec(ipart,nnode),&
            a_b_c%ilen3rec(ipart,nnode),&
            n_cont,inode,ipart,myid,&
            bind_rem,b_rem,lenb_rem,bind,b,&
            a_b_c%istart(ipart,nnode), &
            bmat(1)%mx_abs,parts%mx_mem_grp)
    endif
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
        endif
        jseq=ibseq_rem(nbindbeg+j-1)
        if(jchbeg(jpart)+jseq-1.gt.mx_icover) then
          write(io_lun,912) jchbeg(jpart)+jseq-1,mx_icover
912       format(//'error mat_mult: j_halo arg. out of range',i10,i6)
          icall = 0
          return
        endif
      enddo
    enddo
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
!!    Bug fix for variable function numbers: length now found explicitly from pairs%submat
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
             call cq_abort('mat_tran: itran_addr wrong ', data%itran_addr(ic+nn-1)+count-1,lenat)
          endif
          atrans(data%itran_addr(ic+nn-1):data%itran_addr(ic+nn-1)+count-1)=a_store(posn+1:posn+count)
          posn = posn+count
       enddo
       deallocate(a_store)
    enddo
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
             endif
          enddo ! j=1,nc_nodpart
       endif ! if(nc_nodpart(nn)>0)
    enddo ! nn=1,np_on_node
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
          enddo
          ! Check against A
          if(bmat(nn)%i_part(marker)>amat(nn)%i_part(apos_arr)) then
             cycle ! We've gone past present neighbour - increment A
          else if(bmat(nn)%i_part(marker)==amat(nn)%i_part(apos_arr)) then
             ! Check through seq and part nos
             do while(bmat(nn)%i_seq(marker)<amat(nn)%i_seq(apos_arr).AND.&
                  bmat(nn)%i_part(marker)==amat(nn)%i_part(apos_arr))
                bposn = bposn + nd1*bmat(nn)%ndimj(marker)
                marker = marker+1
             enddo
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
                endif
                marker=marker+1
                bposn = bposn + bsize
                if(bposn-1>lenb) call cq_abort("Length error for B in addscan: ",bposn-1,lenb)
             endif
          endif
          aposn = aposn+bsize
          if(aposn-1>lena) call cq_abort("Length error for A in addscan: ",aposn-1,lena)
       enddo ! nb=1,amat(nn)%n_nab(j)
    endif ! if(amat(nn)%n_nab(j)>0)
  end subroutine addscan
!!***

end module multiply_module
