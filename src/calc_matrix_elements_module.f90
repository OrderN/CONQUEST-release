! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module calc_matrix_elements_module
! ------------------------------------------------------------------------------
! Code area 12: Integration
! ------------------------------------------------------------------------------

!!****h* Conquest/calc_matrix_elements_module *
!!  NAME
!!   calc_matrix_elements_module
!!  PURPOSE
!!   Find matrix elements
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!!   Sometime mid-2000 ?
!!  MODIFICATION HISTORY
!!   17/06/2002 dave
!!    Added headers (ROBODoc for module)
!!   10:38, 06/03/2003 drb 
!!    Corrected array slice pass in act_on_vectors_new (speed problem
!!    under ifc)
!!   10:09, 13/02/2006 drb 
!!    Removed all explicit references to data_ variables and rewrote
!!    in terms of new matrix routines
!!   2006/03/04 06:14 dave
!!    Added use association for function-on-grid data, removed
!!    reference to size
!!   2006/06/15 08:24 dave
!!    Changing in preparation for variable NSF
!!   2008/02/04 16:57 dave
!!    Changed for output to file not stdout
!!   2008/05/23 ast
!!    Added timers
!   2014/09/15 18:30 lat
!!    fixed call start/stop_timer to timer_module (not timer_stdlocks_module !)
!!  SOURCE
!!
module calc_matrix_elements_module

  use global_module,          only: io_lun
  use timer_module,           only: start_timer, stop_timer
  use timer_stdclocks_module, only: tmr_std_integration,     &
                                    tmr_std_allocation,      &
                                    tmr_std_matrices


  implicit none
!!***

contains

!!****f* calc_matrix_elements_module/norb *
!!
!!  NAME 
!!   norb
!!  USAGE
!!   norb(naba_atm structure, atom, block)
!!  PURPOSE
!!   Returns number of orbitals on a given neighbour of a given block
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   T. Miyazaki/D. R. Bowler
!!  CREATION DATE
!!   2006/06/15
!!  MODIFICATION HISTORY
!!   2006/07/05 17:21 dave
!!    Added definition of naba_atm_of_blk
!!  SOURCE
!!
  integer function norb(naba_atm, naba, iprim_blk)

    use naba_blk_module, only: naba_atm_of_blk

    implicit none

    type(naba_atm_of_blk) :: naba_atm
    integer :: naba, iprim_blk

    if(naba == naba_atm%no_of_atom(iprim_blk)) then
       norb = naba_atm%no_of_orb(iprim_blk) - &
              naba_atm%ibeg_orb_atom(naba,iprim_blk) + 1
    else
       norb = naba_atm%ibeg_orb_atom(naba+1,iprim_blk) - &
              naba_atm%ibeg_orb_atom(naba,iprim_blk)
    end if
    
  end function norb
!!***

!!****f* calc_matrix_elements_module/get_matrix_elements_new *
!!
!!  NAME 
!!   get_matrix_elements_new - does integration on grid to get matrix elements
!!  USAGE
!! 
!!  PURPOSE
!!   Loops over blocks on the domain and calculates matrix element partials
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   T. Miyazaki
!!  CREATION DATE
!!   mid-2000
!!  MODIFICATION HISTORY
!!   2006/07/06 08:35 dave
!!    Various preparations and changes for variable NSF
!!   2006/09/29 17:44 dave
!!    Small correction to force size(send_array) into integer
!!   2008/05/23 ast
!!    Added timers
!!   2016/07/20 16:30 nakata
!!    Renamed naba_atm -> naba_atoms_of_blocks
!!  SOURCE
!!
  subroutine get_matrix_elements_new(myid,rem_bucket,matM,gridone,gridtwo)

    use datatypes
    use numbers,             only: zero, one
    use primary_module,      only: domain, bundle
    use naba_blk_module,     only: naba_atm_of_blk
    use bucket_module,       only: local_bucket, remote_bucket
    use GenBlas,             only: gemm
    use comm_array_module,   only: send_array
    use block_module,        only: n_pts_in_block
    use set_blipgrid_module, only: naba_atoms_of_blocks
    use functions_on_grid,   only: gridfunctions, fn_on_grid
    use GenComms,            only: cq_abort
    use memory_module,       only: reg_alloc_mem, reg_dealloc_mem, &
                                   type_dbl
    use global_module,       only: area_integn

    implicit none

    type(remote_bucket), intent(in) :: rem_bucket
    integer,             intent(in) :: myid
    integer,             intent(in) :: gridone,gridtwo
    integer,             intent(in) :: matM ! matM is just ID tag for matrix

    ! Local variables
    real(double), dimension(:), allocatable :: acc_block
    type(local_bucket),    pointer :: loc_bucket
    type(naba_atm_of_blk), pointer :: naba_atm1, naba_atm2

    integer :: iprim_blk, n_dim_one, n_dim_two, i_beg_one, i_beg_two
    integer :: naba1, naba2, ind_halo1, ind_halo2, nsf1, nsf2
    integer :: nonef, ntwof, bucket, where, ii, stat
    integer :: size_send_array

    ! acc_block = zero

    call start_timer(tmr_std_integration)
    ! pointer
    loc_bucket => rem_bucket%locbucket
    naba_atm1  => naba_atoms_of_blocks(loc_bucket%ind_left)
    naba_atm2  => naba_atoms_of_blocks(loc_bucket%ind_right)

    size_send_array = loc_bucket%no_pair_orb
    if(allocated(send_array)) then 
       size_send_array = size(send_array)
       call cq_abort("Error: send_array allocated: ", size_send_array)
    end if
    allocate(send_array(1:size_send_array), STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error allocating send_array in get_matrix_elements: ", &
                       size_send_array)
    call reg_alloc_mem(area_integn, size_send_array, type_dbl)
    send_array(1:size_send_array) = zero

    ! loop over blocks in this node, calculating the contributions 
    ! from each block
    ! we need the following barrier ??
    ! call my_barrier()

    do iprim_blk=1, domain%groups_on_node

       n_dim_one     = naba_atm1%no_of_orb(iprim_blk) !left
       n_dim_two     = naba_atm2%no_of_orb(iprim_blk) !right
       i_beg_one     = (naba_atm1%ibegin_blk_orb(iprim_blk)-1) * n_pts_in_block +1
       i_beg_two     = (naba_atm2%ibegin_blk_orb(iprim_blk)-1) * n_pts_in_block +1
       if(i_beg_one+n_dim_one*n_pts_in_block-1 > gridfunctions(gridone)%size) then
          call cq_abort("get_matrix_elements gridone overflow: ",gridfunctions(gridone)%size, &
               i_beg_one+n_dim_one*n_pts_in_block-1)
       endif
       if(i_beg_two+n_dim_two*n_pts_in_block-1 > gridfunctions(gridtwo)%size) then
          call cq_abort("get_matrix_elements gridtwo overflow: ",gridfunctions(gridtwo)%size, &
               i_beg_two+n_dim_two*n_pts_in_block-1)
       endif

       ! get the overlap contribution from this block

       if(n_dim_two*n_dim_one /= 0) then  ! If the primary block has naba atoms
          allocate(acc_block(n_dim_two*n_dim_one),STAT=stat)
          if(stat/=0) call cq_abort('Error allocating memory in get_matrix_elements: ',n_dim_one,n_dim_two)
          acc_block = zero
          ! for both left and right functions...
          ! To get the contiguous access for acc_block, it should be
          !arranged like acc_block(right,left).
          ! Thus, we use the lower equality (send_array=right*left) in 
          !the following equation.
          !    acc_block(j beta, i alpha)
          !     =\sum_{r_l in block} [ phi_{i alpha} (r_l) * psi_{j beta} (r_l) ]
          !     =\sum_{r_l in block} [ psi_{j beta} (r_l)  * phi_{i alpha} (r_l)]  
          !
          call gemm ( 'T', 'N', n_dim_two, n_dim_one, n_pts_in_block, &
               ONE, gridfunctions(gridtwo)%griddata(i_beg_two:), n_pts_in_block, &
               gridfunctions(gridone)%griddata(i_beg_one:), n_pts_in_block, &
               ZERO, acc_block, n_dim_two )

          ! and accumulate it into the node accumulator

          !---  I HAVE TO BE CAREFUL about WHICH ONE IS LEFT --- 
          !  I have changed the order of naba1 and naba2 to
          ! get contiguous access to send_array.
          ! Note that we can expect loc_bucket%i_h2d(ind_halo2,ind_halo1) changes
          ! gradually in most cases, if ind_halo1 is fixed and both naba2 
          ! and ind_halo2 are orderd in NOPG order
          !
          do naba1=1, naba_atm1%no_of_atom(iprim_blk)    ! left 
             ind_halo1 = naba_atm1%list_atom_by_halo(naba1,iprim_blk)
             do naba2=1, naba_atm2%no_of_atom(iprim_blk) ! right
                ind_halo2 = naba_atm2%list_atom_by_halo(naba2,iprim_blk)
                bucket = loc_bucket%i_h2d(ind_halo2,ind_halo1) !index of the pair
                If(bucket /= 0) then   ! naba1 and naba2 makes pair
                   ii=(bucket-1)
                   ! Note that matrix elements A_{i alpha}{j beta} are arranged 
                   ! like (alpha,beta,ji).  As for ji, j is changing faster than i
                   nonef = norb(naba_atm1,naba1,iprim_blk)
                   ntwof = norb(naba_atm2,naba2,iprim_blk)
                   do nsf2=1, ntwof
                      do nsf1=1, nonef
                         ii = ii+1
                         where = n_dim_two*(naba_atm1%ibeg_orb_atom(naba1, iprim_blk)-1 + nsf1-1) + &
                              naba_atm2%ibeg_orb_atom(naba2, iprim_blk)-1 + nsf2  
                         send_array(ii)=send_array(ii)+acc_block(where)
                      end do
                   end do
                Endif  ! if the two atoms are within a range 

             enddo ! Loop over right naba atoms
          enddo ! Loop over left naba atoms
          deallocate(acc_block)
       endif    ! If the primary block has naba atoms for left & right functions
    end do   ! end loop over primary blocks

    call put_matrix_elements &
         (myid, loc_bucket, rem_bucket, matM )
    deallocate(send_array,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating send_array in get_matrix_elements: ",size_send_array)
    call reg_dealloc_mem(area_integn,size_send_array,type_dbl)
    call stop_timer(tmr_std_integration)
    return
!!***
  contains
    !sbrt put_matrix_elements
    !   30/Jun/2000   Tsuyoshi Miyazaki
    !   Domain responsible node
    !    sends domain partial contribution (send_array).
    !   Bundle responsible node receives and stores them.
    !
    !  Need to be changed if number of functions of each atom varies
    !
    !  20/June/2003 T. MIYAZAKI
    !   We should consider the case where some of the domains have 
    !   no neighbour atoms. (surface with large vacuum area.)
    subroutine put_matrix_elements(myid, loc_bucket, rem_bucket, matM)

      use datatypes
      use numbers
      use mpi
      use bucket_module,     only: local_bucket,remote_bucket
      use dimens,            only: grid_point_volume
      use comm_array_module, only: send_array, recv_array
      use mult_module,       only: matrix_scale, store_matrix_block_pos, &
                                   store_matrix_value_pos, matrix_pos,   &
                                   matrix_index
      use GenComms,          only: cq_abort
      use matrix_module,     only: matrix_halo
      use matrix_data,       only: halo

      implicit none

      type(local_bucket),  intent(in) :: loc_bucket
      type(remote_bucket), intent(in) :: rem_bucket
      integer,             intent(in) :: myid
      integer,             intent(in) :: matM ! matM is just an ID tag for matrix

      real(double), pointer :: recv_ptr(:)
      real(double) :: tmp
      integer :: mynode, myid_ibegin, size_send_array
      integer :: nunit
      integer :: inode, nnd_rem, ibegin, nsize
      integer :: ierr, irc, nsend_req(loc_bucket%mx_recv_node)
      integer :: jnode, ipair, loc1, loc2, ii, isf1, isf2, stat, msize
      integer, dimension(MPI_STATUS_SIZE) :: nrecv_stat, nwait_stat
      integer, parameter :: tag = 1
            
      mynode = myid + 1
      ierr = 0 
      size_send_array = size(send_array)
      call matrix_scale(zero, matM)
      if (loc_bucket%no_recv_node > 0) then
         do inode = 1, loc_bucket%no_recv_node ! Loop over recv nodes
            nnd_rem = loc_bucket%list_recv_node(inode)
            if (inode == 1) then
               ibegin = 1
            else
               ibegin = loc_bucket%ibegin_pair_orb(inode)
            end if
            if (inode == loc_bucket%no_recv_node) then
               nsize = loc_bucket%no_pair_orb - &
                       loc_bucket%ibegin_pair_orb(inode) + 1
            else
               nsize = loc_bucket%ibegin_pair_orb(inode + 1) - &
                       loc_bucket%ibegin_pair_orb(inode)
            end if
            !SENDING or KEEPING the First Address 
            if (nnd_rem == mynode) then
               myid_ibegin = ibegin
            else if (nsize > 0) then
               call MPI_issend(send_array(ibegin), nsize,            &
                               MPI_DOUBLE_PRECISION, nnd_rem-1, tag, &
                               MPI_COMM_WORLD, nsend_req(inode), ierr)
               if (ierr /= 0) &
                    call cq_abort('ERROR in MPI_issend in put_matrix_elements', &
                                  ierr)
            end if
         end do
      end if ! (loc_bucket%no_recv_node >0) 

      if (rem_bucket%no_send_node > 0) then
         msize = 0
         do jnode=1,rem_bucket%no_send_node  ! Loop over sending nodes
            msize = max(msize,rem_bucket%no_of_pair_orbs(jnode))
         end do
         allocate(recv_array(msize),STAT=stat)
         if (stat/=0) call cq_abort("Error allocating recv_array in put_matrix_elements: ",msize)
         call reg_alloc_mem(area_integn,msize,type_dbl)
         call start_timer(tmr_std_matrices)
         do jnode=1,rem_bucket%no_send_node  ! Loop over sending nodes
            nnd_rem=rem_bucket%list_send_node(jnode)
            nsize=rem_bucket%no_of_pair_orbs(jnode)
            if (nsize>0) then
               recv_array = zero
               if (nnd_rem == mynode) then
                  if (myid_ibegin+nsize-1<=size_send_array) then
                     recv_ptr => send_array(myid_ibegin:myid_ibegin+nsize-1)
                  else
                     call cq_abort("Overflow in put_matrix_elements: ",myid_ibegin+nsize-1,size_send_array)
                  end if
               else
                  call MPI_recv(recv_array,nsize,MPI_DOUBLE_PRECISION,nnd_rem-1, &
                       tag,MPI_COMM_WORLD,nrecv_stat,ierr)
                  if (ierr /= 0) call cq_abort('ERROR in MPI_recv in put_matrix_elements',ierr)
                  recv_ptr => recv_array
               end if
               do ipair=1,rem_bucket%no_of_pair(jnode)
                  loc2=rem_bucket%bucket(ipair,jnode)%ibegin_pair!(ipair-1)*nunit
                  loc1=matrix_pos(matM,rem_bucket%bucket(ipair,jnode)%iprim,rem_bucket%bucket(ipair,jnode)%jhalo)
                  if (loc1 /= 0) then
                     !call store_matrix_block_pos(matM,loc1,recv_ptr(loc2+1:loc2+nunit)*grid_point_volume,nunit)
                     ! VarNSF: change nunit to use orbitals from iprim/jhalo
                     ii=0
                     do isf2 = 1,halo(matrix_index(matM))%ndimj(rem_bucket%bucket(ipair,jnode)%jhalo)
                        do isf1 = 1,halo(matrix_index(matM))%ndimi(rem_bucket%bucket(ipair,jnode)%iprim)
                           tmp = recv_ptr(loc2+ii)*grid_point_volume
                           call store_matrix_value_pos(matM,loc1+ii,tmp)
                           ii=ii+1
                        end do
                     end do
                  else 
                     write(io_lun,*) 'ERROR? in remote_bucket',ipair,myid
                     write(io_lun,*) ' a pair of two atoms in a domain partial node &
                          & is not a neighbour pair in a bundle node. ??? '
                  end if
               end do        ! Loop over sent pairs
            end if
         end do        ! Loop over sending nodes
         call stop_timer(tmr_std_matrices)
         deallocate(recv_array,STAT=stat)
         if (stat/=0) call cq_abort("Error deallocating recv_array in put_matrix_elements: ",msize)
         call reg_dealloc_mem(area_integn,msize,type_dbl)
      end if !(rem_bucket%no_send_node >0) 

      if (loc_bucket%no_recv_node >0) then
         do inode=1,loc_bucket%no_recv_node
            if (inode == loc_bucket%no_recv_node) then
               nsize=loc_bucket%no_pair_orb-loc_bucket%ibegin_pair_orb(inode)+1
            else
               nsize=loc_bucket%ibegin_pair_orb(inode+1) - loc_bucket%ibegin_pair_orb(inode)
            end if
            nnd_rem=loc_bucket%list_recv_node(inode)
            if (nnd_rem /= mynode.AND.nsize>0) then
               call MPI_WAIT(nsend_req(inode),nwait_stat,ierr)
               if (ierr /= 0) call cq_abort('ERROR in MPI_WAIT in put_matrix_elements',ierr)
            end if
         end do
      end if  ! (loc_bucket%no_recv_node >0)
      return

    end subroutine put_matrix_elements

  end subroutine get_matrix_elements_new

!--------------------------------------------------------------------------
!sbrt act_on_vectors_new&
!    (myid,rem_bucket,data_M,gridone,gridtwo)
!------------ NEW VERSION -------------------------------------------------
! Now we assume, every atom has same NONEF(NTWOF) for 
! gridone(gridtwo).
! We have to change 
!     nonef -> nonef(ia) ia: naba atm
!     calculation of ibegin_blk in set_blipgrid_module
!
! The subroutine accumulates for gridone
! it calculates gridone = gridone + matM gridtwo
!
! Modifications:
!  2016/07/20 16:30 nakata
!   Renamed naba_atm -> naba_atoms_of_blocks
  subroutine act_on_vectors_new(myid,rem_bucket,matM,gridone,gridtwo)
    ! Modules and Dummy arguments
    use datatypes
    use numbers
    use primary_module,    only:domain
    use naba_blk_module,   only:naba_atm_of_blk
    use set_blipgrid_module, only: naba_atoms_of_blocks
    use bucket_module,     only:local_bucket,remote_bucket
    use GenBlas,           only:axpy
    use comm_array_module, only:send_array
    use block_module,      only:n_pts_in_block
    use functions_on_grid, only: gridfunctions, fn_on_grid
    use GenComms, only: cq_abort
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use global_module, only: area_integn

    implicit none

    !type(naba_atm_of_blk),intent(in),target:: naba_atm1,naba_atm2
    !type(local_bucket),intent(in)          :: loc_bucket
    type(remote_bucket), target, intent(in) :: rem_bucket
    integer,                     intent(in) :: myid
    integer :: matM
    integer :: gridone, gridtwo   ! a set of left and right functions

    ! Local variables
    type(naba_atm_of_blk), pointer :: naba_atm1, naba_atm2
    type(local_bucket),    pointer :: loc_bucket
    integer      :: iprim_blk, n_dim_one, n_dim_two
    integer      :: naba1, naba2, ind_halo1, ind_halo2, bucket
    integer      :: nonef, ntwof, ind1,ind2
    integer      :: nsf1, nsf2, ii, stat
    real(double) :: factor_M

    call start_timer(tmr_std_integration)
    !(pointer)
    !  This structure is introduced to keep the strict connection
    !  of remote_bucket <-> local_bucket <-> naba_atm_of_blk.
    loc_bucket => rem_bucket%locbucket
    naba_atm1  => naba_atoms_of_blocks(loc_bucket%ind_left)
    naba_atm2  => naba_atoms_of_blocks(loc_bucket%ind_right)

    allocate(send_array(loc_bucket%no_pair_orb), STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error allocating send_array in act_on_vectors: ", &
                       loc_bucket%no_pair_orb)
    call reg_alloc_mem(area_integn, loc_bucket%no_pair_orb, type_dbl)
    ! collects matrix elements and store them in send_array
    call collect_matrix_elements(myid, loc_bucket, rem_bucket, matM)

    ! Initialization of output (gridone)
    !      gridone(:) = zero     ! NOT INITIALISE
    ! FOR get_non_local_gradient !
    !  19/Sep/00 Tsuyoshi Miyazaki
    ! loop over blocks in this node
    do iprim_blk=1, domain%groups_on_node

       !  In the future we have to prepare n_dim_one & n_dim_two 
       ! for each iprim_blk by counting all NONEF and NTWOF for 
       ! the neighbour atoms of the block.
       !  In this case, the loops ovre nsf1 & nsf2 have to be
       ! deleted.
       n_dim_one = naba_atm1%no_of_orb(iprim_blk) !left
       n_dim_two = naba_atm2%no_of_orb(iprim_blk) !right

       ! IN (get_matrix_elements)
       !  if (n_dim_two.ne.0)  &
       !   call gemm (...)

       if (n_dim_one*n_dim_two /= 0) then
          ! I HAVE TO BE CAREFUL about WHICH ONE IS LEFT
          do naba1=1, naba_atm1%no_of_atom(iprim_blk)    ! left atom
             ind_halo1 = naba_atm1%list_atom_by_halo(naba1,iprim_blk)
             if(ind_halo1 > loc_bucket%no_halo_atom1) &
                  call cq_abort('ERROR in no_of_halo_atoms for left',ind_halo1,loc_bucket%no_halo_atom1)
             do naba2=1, naba_atm2%no_of_atom(iprim_blk) ! right atom
                ind_halo2 = naba_atm2%list_atom_by_halo(naba2,iprim_blk)
                if(ind_halo2 > loc_bucket%no_halo_atom2) &
                     call cq_abort('ERROR in no_of_halo_atoms for right',ind_halo2,loc_bucket%no_halo_atom2)

                bucket = loc_bucket%i_h2d(ind_halo2,ind_halo1) !index of the pair
                if(bucket > loc_bucket%no_pair_orb) &
                     call cq_abort('ERROR : bucket in act_on_vectors_new',bucket,loc_bucket%no_pair_orb)
                If(bucket /= 0) then   ! naba1 and naba2 makes pair
                   ii=bucket-1
                   nonef = norb(naba_atm1, naba1, iprim_blk)
                   ntwof = norb(naba_atm2, naba2, iprim_blk)
                   do nsf2=1, ntwof
                      ind2=n_pts_in_block*(naba_atm2%ibegin_blk_orb(iprim_blk)-1+ &
                           naba_atm2%ibeg_orb_atom(naba2,iprim_blk)-1+(nsf2-1))+1
                      do nsf1=1, nonef                   
                         ind1=n_pts_in_block*(naba_atm1%ibegin_blk_orb(iprim_blk)-1+ &
                              naba_atm1%ibeg_orb_atom(naba1,iprim_blk)-1+(nsf1-1))+1
                         ii = ii+1
                         factor_M=send_array(ii)
                         call axpy(n_pts_in_block, factor_M, &
                              gridfunctions(gridtwo)%griddata(ind2:ind2+n_pts_in_block-1), 1, &
                              gridfunctions(gridone)%griddata(ind1:ind1+n_pts_in_block-1), 1)
                      end do
                   end do
                Endif  ! if the two atoms are within a range 

             enddo ! Loop over right naba atoms
          enddo ! Loop over left naba atoms
       endif  ! If the primary block has naba atoms
    end do   ! end loop over primary blocks
    deallocate(send_array,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating send_array in act_on_vectors: ",loc_bucket%no_pair_orb)
    call reg_dealloc_mem(area_integn,loc_bucket%no_pair_orb,type_dbl)
    call stop_timer(tmr_std_integration)

    return
  contains
    !-----------------------------------------------------------------
    !sbrt collect_matrix_elements
    !  Created : 3/Jul/2000 
    !             by Tsuyoshi Miyazaki
    !
    !    This subroutine is called by (act_on_vector_new).
    !    Bundle responsible nodes send their responsible matrix elements
    !   to the domain responsible nodes, which will calculate
    !   \sum_{j} [A_{i alpha; j beta} phi_{j beta} (r_l)].
    !
    !    I am sorry if this subroutine irritates you to read.
    !   I recommend you to read (put_matrix_elements) before reading 
    !   this subroutine.
    !    This subroutine will do inverse things of (put_matrix_elements),
    !   which is called by (get_matrix_elements_new).
    !   <send_array>, which will be used to store their receiving arrays
    !   in this subroutine, has the same size as <send_array> 
    !   in (put_matrix_elements). 
    !
    !    Need to be changed if number of functions of each atom varies
    !
    !  20/June/2003 T. MIYAZAKI
    !   We should consider the case where some of the domains have 
    !   no neighbour atoms. (surface with large vacuum area.)
    subroutine collect_matrix_elements(myid, loc_bucket, rem_bucket, matM)

      use datatypes
      use mpi
      use bucket_module,     only: local_bucket,remote_bucket
      use comm_array_module, only: send_array, recv_array
      use GenComms,          only: my_barrier, cq_abort
      use mult_module,       only: return_matrix_value_pos, &
                                   matrix_pos, matrix_index
      use matrix_data,       only: halo

      implicit none

      type(local_bucket), intent(in) :: loc_bucket
      type(remote_bucket),intent(in) :: rem_bucket
      integer, intent(in) :: myid
      integer :: matM

      real(double),pointer :: recv_ptr(:)
      integer :: mynode,myid_ibegin
      integer :: inode,nnd_rem,ibegin,nsize,msize
      integer :: jnode,ipair,loc1,loc2,ii,isf1,isf2, off
      integer :: ierr,irc,stat
      integer :: nrecv_req(loc_bucket%mx_recv_node)
      integer :: nwait_stat(MPI_STATUS_SIZE)
      integer,parameter :: tag=1

      mynode=myid+1
      ierr=0 

      ! First prepare receiving the array AND keep the first address
      ! of the data within my node in send_array.
      ! 
      if(loc_bucket%no_recv_node > 0) then
         do inode=1,loc_bucket%no_recv_node ! Loop over recv nodes
            nnd_rem=loc_bucket%list_recv_node(inode)
            if(inode ==1) then
               ibegin=1
            else
               ibegin=loc_bucket%ibegin_pair_orb(inode)
            endif
            if(inode == loc_bucket%no_recv_node) then
               nsize=loc_bucket%no_pair_orb-loc_bucket%ibegin_pair_orb(inode)+1
            else
               nsize=loc_bucket%ibegin_pair_orb(inode+1) - loc_bucket%ibegin_pair_orb(inode)
            endif

            ! Prepare to receive OR keep the First Address
            if(nnd_rem == mynode) then
               myid_ibegin=ibegin
            else if(nsize>0) then
               call MPI_irecv(send_array(ibegin),nsize,MPI_DOUBLE_PRECISION,&
                    nnd_rem-1,tag,MPI_COMM_WORLD,nrecv_req(inode),ierr)
               if(ierr /= 0) call cq_abort('ERROR in MPI_irecv in collect_matrix_elements',ierr)
            endif
         enddo          ! Loop over receiving nodes
      endif  ! (loc_bucket%no_recv_node > 0) 

      if(rem_bucket%no_send_node > 0) then
         msize=0
         do jnode=1,rem_bucket%no_send_node  ! Loop over receiving nodes
            msize=max(msize,rem_bucket%no_of_pair_orbs(jnode))
         end do
         allocate(recv_array(msize),STAT=stat)
         if(stat/=0) call cq_abort("Error allocating recv_array in collect_matrix_elements: ",msize)
         call reg_alloc_mem(area_integn,msize,type_dbl)
         off = 0
         call start_timer(tmr_std_matrices)
         do jnode=1,rem_bucket%no_send_node  ! Loop over receiving nodes
            recv_array = zero
            nnd_rem=rem_bucket%list_send_node(jnode)
            nsize=rem_bucket%no_of_pair_orbs(jnode)
            if(nsize>0) then
               if(nnd_rem == mynode) then
                  recv_ptr => send_array(myid_ibegin:myid_ibegin+nsize-1)
               else
                  recv_ptr => recv_array
               endif
               !
               !Making recv_ptr (actually, this is the array which will be sent)
               !
               do ipair=1,rem_bucket%no_of_pair(jnode)
                  loc2=rem_bucket%bucket(ipair,jnode)%ibegin_pair
                  loc1=matrix_pos(matM,rem_bucket%bucket(ipair,jnode)%iprim,rem_bucket%bucket(ipair,jnode)%jhalo)
                  if(loc1 /= 0) then 
                     ii=0
                     do isf2 = 1,halo(matrix_index(matM))%ndimj(rem_bucket%bucket(ipair,jnode)%jhalo)
                        do isf1 = 1,halo(matrix_index(matM))%ndimi(rem_bucket%bucket(ipair,jnode)%iprim)
                           recv_ptr(loc2+ii)= return_matrix_value_pos(matM,loc1+ii)
                           ii=ii+1
                        end do
                     end do
                  else
                     write(io_lun,*) 'ERROR? of remote_bucket in &
                          & collect_matrix_elements',ipair,myid
                     write(io_lun,*) ' a pair of two atoms in a domain partial node &
                          & is not a neighbour pair in a bundle node. ??? '
                  endif
               enddo        ! Loop over sending pairs

               if(nnd_rem /= mynode) then
                  call MPI_ssend(recv_array,nsize,MPI_DOUBLE_PRECISION,nnd_rem-1, &
                       tag,MPI_COMM_WORLD,ierr)
                  if(ierr /= 0) call cq_abort('ERROR in MPI_ssend in put_matrix_elements',ierr)
               endif
               !deallocate(recv_array,STAT=stat)
               !if(stat/=0) call cq_abort("Error deallocating recv_array in collect_matrix_elements: ",nsize)
               off = off+nsize
            end if
         enddo        ! Loop over sending nodes
         call stop_timer(tmr_std_matrices)
      endif   !(rem_bucket%no_send_node > 0) 

      !Check whether MPI_irecv has been finished.
      if(loc_bucket%no_recv_node > 0) then
         do inode=1,loc_bucket%no_recv_node
            if(inode == loc_bucket%no_recv_node) then
               nsize=loc_bucket%no_pair_orb-loc_bucket%ibegin_pair_orb(inode)+1
            else
               nsize=loc_bucket%ibegin_pair_orb(inode+1) - loc_bucket%ibegin_pair_orb(inode)
            endif
            nnd_rem=loc_bucket%list_recv_node(inode)
            if(nnd_rem /= mynode.AND.nsize>0) then
               call MPI_WAIT(nrecv_req(inode),nwait_stat,ierr)
               if(ierr /= 0) call cq_abort('ERROR in MPI_WAIT in put_matrix_elements',ierr)
            endif
         enddo
      endif   ! (loc_bucket%no_recv_node > 0)
      if(rem_bucket%no_send_node>0) then
         if(allocated(recv_array)) then
            nullify(recv_ptr)
            deallocate(recv_array,STAT=stat)
            if(stat/=0) call cq_abort("Error deallocating recv_array in collect_matrix_elements: ",nsize)
            call reg_dealloc_mem(area_integn,msize,type_dbl)
         else
            write(io_lun,*) 'Wierd: recv_array no allocated'
         endif
      endif
      return
    end subroutine collect_matrix_elements
  end subroutine act_on_vectors_new
end module calc_matrix_elements_module
