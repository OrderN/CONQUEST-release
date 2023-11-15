! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module block_module
! ------------------------------------------------------------------------------
! Code area 8: indexing
! ------------------------------------------------------------------------------

!!****h* Conquest/block_module *
!!  NAME
!!   block_module
!!  PURPOSE
!!     This module shows the numbers
!!    of integration grid points in a block along
!!    x, y and z direction.
!!  (In the old version, they are in_block_x,
!!     in_block_y, and in_block_z)
!!
!!     In the future, this module should be merged into 
!!    group_module or primary_module, because integration
!!    grid points are just the members of a primary set,
!!    domain.
!!  
!!     It should be also noted that the order of integration
!!    grid points in a block is assumed to be the following,
!!            order=0
!!       do nz=1,nz_in_block
!!        do ny=1,ny_in_block
!!         do nx=1,nx_in_block
!!            order=order+1
!!         enddo
!!        enddo
!!       enddo
!! 
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!!   11/07/2000
!!  MODIFICATION HISTORY
!!   15/03/2002 dave
!!    Added ROBODoc header, RCS id tags
!!   18/03/2002 dave
!!    Trivial tidying
!!   2004/11/10 drb
!!    Added set_blocks_from_new (to take block structure and create old-style data)
!!   2006/10/10 16:49 dave
!!    Removed grid_point_volume (to be found in dimens.module.f90)
!!   2008/02/04 08:30 dave
!!    Changed for output to file not stdout
!!   2008/05/15 ast
!!    Added some timers
!!   2014/09/15 18:30 lat
!!    fixed call start/stop_timer to timer_module (not timer_stdlocks_module !)
!!  SOURCE
!!
module block_module

  use datatypes
  use global_module,          only: io_lun
  use timer_module,           only: start_timer,stop_timer
  use timer_stdclocks_module, only: tmr_std_indexing,tmr_std_allocation

  implicit none
  save

  integer :: nx_in_block,ny_in_block,nz_in_block
  integer :: n_pts_in_block
  integer :: n_blocks
  integer :: in_block_x,in_block_y,in_block_z
  !real(double) :: grid_point_volume

  integer, parameter :: blocks_file    = 0
  integer, parameter :: blocks_raster  = 1
  integer, parameter :: blocks_hilbert = 2

!!***

contains

!!****f* block_module/set_blocks_from_new *
!!
!!  NAME 
!!   set_blocks_from_new
!!  USAGE
!! 
!!  PURPOSE
!!   Sets up the old-style block information from a new-style initialisation (temporary
!!   measure until I can rewrite the FFT routines !)
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   2004/11/10
!!  MODIFICATION HISTORY
!!   2008/05/16 ast
!!    Added timers
!!  SOURCE
!!
  subroutine set_blocks_from_new()

  use grid_index, ONLY: ind_block_x, ind_block_y, ind_block_z
  use GenComms, ONLY: inode, cq_abort
  use group_module, ONLY: blocks
  use maxima_module, ONLY: maxblocks

  ! Local variables
  integer :: icount, iblock, ind_group, stat

  call start_timer(tmr_std_indexing)

  call start_timer(tmr_std_allocation)
  allocate(ind_block_x(maxblocks),ind_block_y(maxblocks),ind_block_z(maxblocks),STAT=stat)
  call stop_timer(tmr_std_allocation)
  if(stat/=0) call cq_abort("Error allocating ind_block_x: ",maxblocks,stat)
  icount = 0
  n_blocks = blocks%ng_on_node(inode)
  do iblock = 1,blocks%ng_on_node(inode)
     icount = icount + 1
     ind_group = blocks%ngnode(blocks%inode_beg(inode)+icount-1)
     ind_block_x(icount) = 1+(ind_group-1)/(blocks%ngcelly*blocks%ngcellz)
     ind_block_y(icount) = 1+(ind_group-1-(ind_block_x(icount)-1)*blocks%ngcelly*blocks%ngcellz)/blocks%ngcellz
     ind_block_z(icount) = ind_group - &
          (ind_block_x(icount)-1)*blocks%ngcelly*blocks%ngcellz - (ind_block_y(icount)-1)*blocks%ngcellz
  end do
  call stop_timer(tmr_std_indexing)
  end subroutine set_blocks_from_new
!!***

!!****f* block_module/set_blocks_from_old *
!!
!!  NAME 
!!   set_blocks_from_old
!!  USAGE
!! 
!!  PURPOSE
!!   Nakes the type(group_set) :: blocks
!!
!!   Rasters over grid blocks and distributes them to nodes: default behaviour, though not 
!!   necessarily the best
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   T. Miyazaki
!!  CREATION DATE
!!   12/07/2000
!!  MODIFICATION HISTORY
!!   10:08, 29/10/2002 dave
!!    Added ROBODoc headers
!!   11:25, 12/11/2004 dave
!!    Removed references to distribute_atom_module
!!   2006/08/03 08:04 dave
!!    Incorporated into block_module
!!   2006/10/10 16:59 dave
!!    in_block_ variables now part of module, so removed from argument list
!!   2007/04/18 17:26 dave
!!    Added block assignment following partitions and choice flag
!!   2007/05/18 15:17 dave
!!    Added lines to write out block distribution for Hilbert partitioner
!!   2008/01/21 Veronika
!!    Changed write format "Minimum blocs/proc ..."
!!   2008/05/16 ast
!!    Added timers
!!   2008/09/01 08:21 dave
!!    Added io_ routines from input_module
!!  SOURCE
!!
  subroutine set_blocks_from_old(n_grid_x, n_grid_y, n_grid_z)

    !Modules and Dummy Arguments
    use datatypes
    use numbers
    use global_module, ONLY: numprocs, iprint_index, flag_assign_blocks
    use construct_module, ONLY: init_group
    use group_module,  ONLY: blocks, allocate_group_set, make_cc2, parts
    use maxima_module, ONLY: maxngrid, maxblocks
    use GenComms, ONLY: inode, ionode, cq_abort
    use input_module, ONLY: io_assign, io_close

    implicit none

    integer,intent(in):: n_grid_x, n_grid_y, n_grid_z

    !Local variables
    integer:: n_block_x, n_block_y, n_block_z,   &
         n_per_proc, n_first, n_all_blocks,              &
         nfb_x, nfb_y, nfb_z, nc_x, nc_y, nc_z,  &
         n_add_x, n_add_y, n_add_z, nx, ny, nz, icount, lun


    integer:: mx_gedge_tmp,ind_block,nnd, maxtmp, i, cx, cy, cz, cc
    integer:: minblocks

    integer, dimension(:), allocatable :: blocks_per_proc, proc_block, proc_which_block

    call start_timer(tmr_std_indexing)
    !--- find numbers of blocks in each direction
    n_block_x = n_grid_x/in_block_x
    n_block_y = n_grid_y/in_block_y
    n_block_z = n_grid_z/in_block_z
    if(n_block_x*in_block_x/=n_grid_x) call cq_abort("Block and grid sizes incompatible ! ",in_block_x,n_grid_x)
    if(n_block_y*in_block_y/=n_grid_y) call cq_abort("Block and grid sizes incompatible ! ",in_block_y,n_grid_y)
    if(n_block_z*in_block_z/=n_grid_z) call cq_abort("Block and grid sizes incompatible ! ",in_block_z,n_grid_z)
    blocks%ngcellx=n_block_x
    blocks%ngcelly=n_block_y
    blocks%ngcellz=n_block_z
    !-- set variables in block_module
    nx_in_block=in_block_x
    ny_in_block=in_block_y
    nz_in_block=in_block_z
    mx_gedge_tmp=max(n_block_x,n_block_y,n_block_z)
    if(flag_assign_blocks==blocks_raster) then
       !--- distribution of blocks over nodes:
       !--- the first  n_first  nodes each have  n_per_proc + 1  blocks;
       !--- the remaining  numprocs - n_first  nodes each have  n_per_proc  blocks.
       n_all_blocks = n_block_x * n_block_y * n_block_z
       n_per_proc = n_all_blocks / numprocs
       if(n_per_proc==0) call cq_abort("There must be at least one block per node")
       n_first = n_all_blocks - n_per_proc * numprocs
       if(inode==ionode.AND.iprint_index>0) then
          if (n_first.ne.0) then
             write(io_lun,fmt='(10x,"The first ",i5," processors each have ",i5," blocks"/&
                  &10x,"The remaining ",i5," processors each have ",i5," blocks.")') &
                  n_first,n_per_proc+1,numprocs-n_first,n_per_proc
          else
             write(io_lun,fmt='(10x,"All processors have ",i5," blocks.")') n_per_proc
          end if
       endif
       ! Define maximum number of grid points on any processor for FFTs
       maxblocks = n_per_proc+1
       maxngrid = n_pts_in_block * (n_per_proc+1)
       ! FFT check
       maxtmp = n_grid_x * (n_grid_y*n_grid_z/numprocs + 1)
       if(maxtmp>maxngrid) maxngrid = maxtmp
       maxtmp = n_grid_y * (n_grid_z*n_grid_x/numprocs + 1)
       if(maxtmp>maxngrid) maxngrid = maxtmp
       maxtmp = n_grid_z * (n_grid_x*n_grid_y/numprocs + 1)
       if(maxtmp>maxngrid) maxngrid = maxtmp
       call init_group(blocks, maxblocks, mx_gedge_tmp, n_all_blocks, n_pts_in_block, numprocs)

       !--- number of blocks on current node
       blocks%ng_on_node(:)=0
       do nnd=1, numprocs
          if(nnd .gt. n_first) then
             blocks%ng_on_node(nnd) = n_per_proc
          else
             blocks%ng_on_node(nnd) = n_per_proc + 1
          endif
       enddo
       !--- index of first block on current node
       blocks%inode_beg=0
       do nnd=1, numprocs
          if(nnd == 1) then
             blocks%inode_beg(nnd)=1
          else
             blocks%inode_beg(nnd)= &
                  blocks%inode_beg(nnd-1)+blocks%ng_on_node(nnd-1)
          endif
       enddo
       !--- make list of indices of blocks on current node
       blocks%ngnode(:)=0
       blocks%inv_ngnode(:)=0

       nfb_x = 1 ; nfb_y = 1 ; nfb_z = 1
       nc_x = 1  ; nc_y = 1  ; nc_z = 1
       icount=0

       n_add_z = 0
       do nz = 1, n_block_z
          nc_z = nc_z + n_add_z
          n_add_z = nfb_z

          n_add_y = 0
          do ny = 1, n_block_y
             nc_y = nc_y + n_add_y
             n_add_y = nfb_y

             n_add_x = 0
             do nx = 1, n_block_x
                nc_x = nc_x + n_add_x
                n_add_x = nfb_x
                icount=icount+1
                ind_block = nc_z+(nc_y-1)*n_block_z &
                     +(nc_x-1)*n_block_y*n_block_z

                blocks%ngnode(icount) = ind_block
                blocks%inv_ngnode(ind_block)=icount
             enddo
             nfb_x = -nfb_x

          enddo
          nfb_y = -nfb_y

       enddo
       if(inode==ionode) then
          call io_assign(lun)
          open(unit=lun,file='raster_make_blk.dat')
          write(lun,fmt='(3i8)') n_block_x, n_block_y, n_block_z
          write(lun,fmt='(i8)') numprocs
          do nnd = 1,numprocs
             write(lun,fmt='(3i8)') nnd,blocks%ng_on_node(nnd),blocks%inode_beg(nnd)
             do nx = 1,blocks%ng_on_node(nnd)
                write(lun,fmt='(2i8)') nx,blocks%ngnode(blocks%inode_beg(nnd)+nx-1)
             end do
          end do
          call io_close(lun)
       end if
    else if(flag_assign_blocks==blocks_hilbert) then
       n_all_blocks = n_block_x * n_block_y * n_block_z
       allocate(blocks_per_proc(numprocs),proc_block(n_all_blocks),proc_which_block(n_all_blocks))
       ! This counts blocks per processor
       blocks_per_proc = 0
       ! This stores the processor for each block
       proc_block = 0
       ! This stores WHICH block on the processor each block is
       proc_which_block = 0
       maxblocks = 0
       i=0
       ! Work out max blocks on a processor: loop over blocks
       do nz = 1, n_block_z
          do ny = 1, n_block_y
             do nx = 1, n_block_x
                i = i+1
                ! Convert centre of block into a partition index
                cx = floor(parts%ngcellx*(real(nx,double)-0.5)/real(n_block_x,double))
                cy = floor(parts%ngcelly*(real(ny,double)-0.5)/real(n_block_y,double))
                cz = floor(parts%ngcellz*(real(nz,double)-0.5)/real(n_block_z,double))
                cc = 1+cz+parts%ngcellz*(cy+parts%ngcelly*cx)
                proc_block(i) = parts%i_cc2node(cc)
                blocks_per_proc(proc_block(i)) = blocks_per_proc(proc_block(i))+1
                if(blocks_per_proc(proc_block(i))>maxblocks) maxblocks = blocks_per_proc(proc_block(i))
                proc_which_block(i) = blocks_per_proc(proc_block(i))
                !write(io_lun,*) 'Block: ',i,proc_block(i),proc_which_block(i)
             end do
          end do
       end do
       ! Define maximum number of grid points on any processor for FFTs
       maxngrid = n_pts_in_block * maxblocks
       ! FFT check
       maxtmp = n_grid_x * (n_grid_y*n_grid_z/numprocs + 1)
       if(maxtmp>maxngrid) maxngrid = maxtmp
       maxtmp = n_grid_y * (n_grid_z*n_grid_x/numprocs + 1)
       if(maxtmp>maxngrid) maxngrid = maxtmp
       maxtmp = n_grid_z * (n_grid_x*n_grid_y/numprocs + 1)
       if(maxtmp>maxngrid) maxngrid = maxtmp
       call init_group(blocks, maxblocks, mx_gedge_tmp, n_all_blocks, n_pts_in_block, numprocs)
       !--- number of blocks on current node
       blocks%ng_on_node(:)=0
       minblocks = 1e8
       do nnd=1, numprocs
          blocks%ng_on_node(nnd) = blocks_per_proc(nnd)
          if(blocks_per_proc(nnd)<minblocks) minblocks = blocks_per_proc(nnd)
       enddo
       if(inode==ionode.AND.iprint_index>1) then
          write(io_lun,fmt='(4x,"Minimum blocs/proc is ",i8,". Maximum blocs/proc is ",i8)') minblocks, maxblocks
       end if
       !--- index of first block on current node
       blocks%inode_beg=0
       do nnd=1, numprocs
          if(nnd == 1) then
             blocks%inode_beg(nnd)=1
          else
             blocks%inode_beg(nnd)= &
                  blocks%inode_beg(nnd-1)+blocks%ng_on_node(nnd-1)
          endif
       enddo
       ! ngnode gives the CC label of block as we loop over them in NODE ORDER
       ! Here we loop over blocks in CC order and work out the node order label
       i=0
       do nz = 1, n_block_z
          do ny = 1, n_block_y
             do nx = 1, n_block_x
                i = i+1
                nnd = proc_block(i)
                icount = blocks%inode_beg(nnd) + proc_which_block(i)-1
                ind_block = nz+blocks%ngcellz*(ny-1+blocks%ngcelly*(nx-1))
                blocks%ngnode(icount) = ind_block
                blocks%inv_ngnode(ind_block)=icount
             end do
          end do
       end do
       deallocate(blocks_per_proc,proc_block,proc_which_block)
       if(inode==ionode) then
          call io_assign(lun)
          open(unit=lun,file='hilbert_make_blk.dat')
          write(lun,fmt='(3i8)') n_block_x, n_block_y, n_block_z
          write(lun,fmt='(i8)') numprocs
          do nnd = 1,numprocs
             write(lun,fmt='(3i8)') nnd,blocks%ng_on_node(nnd),blocks%inode_beg(nnd)
             do nx = 1,blocks%ng_on_node(nnd)
                write(lun,fmt='(2i8)') nx,blocks%ngnode(blocks%inode_beg(nnd)+nx-1)
             end do
          end do
          call io_close(lun)
       end if
    else 
       call cq_abort("Unknown block-assignment mode: ",flag_assign_blocks)
    end if
    !--- member information     -- SHOULD BE MEANINGLESS
    blocks%nm_group(:)=0
    blocks%icell_beg(:)=0
    do ind_block=1,n_all_blocks
       blocks%nm_group(ind_block)=n_pts_in_block
       blocks%icell_beg(ind_block)=(ind_block-1)*n_pts_in_block+1
    enddo
    call stop_timer(tmr_std_indexing)
    !--- make_cc2
    ! This already times for indexing
    ! That's why the timer is stopped above
    call make_cc2(blocks,numprocs)

    return
  end subroutine set_blocks_from_old
!!***

!!****f* block_module/set_domains *
!!
!!  NAME 
!!   set_domains
!!  USAGE
!! 
!!  PURPOSE
!!   sets the necessary arrays to define the domain
!!   corresponding to each PE. Present version
!!   is based on organisation of grid-points into
!!   blocks. Input for this routine is generated by
!!   a prior call of the subroutine blocks.
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   M.J.Gillan
!!  CREATION DATE
!!   19/03/96
!!  MODIFICATION HISTORY
!!   25/05/2001 dave
!!    Stripped subroutine call, added ROBODoc header
!!   20/06/2001 dave
!!    Added RCS Id and Log tags and cq_abort
!!   2006/11/21 17:20 dave
!!    Included in block_module
!!   2008/05/16 ast
!!    Added timer
!!  SOURCE
!!
  subroutine set_domains

    !use dimens, ONLY: n_my_grid_points
    use grid_index, ONLY: ind_block_x, ind_block_y, ind_block_z, &
         grid_point_x, grid_point_y, grid_point_z, &
         grid_point_block, grid_point_position
    use GenComms, ONLY: cq_abort

    implicit none

    ! Shared variables

    ! Local variables
    integer :: n_pts, nb,  &
         n0_x, n0_y, n0_z, nx, ny, nz

    call start_timer(tmr_std_indexing)
    ! This subroutine replaces a previous version in which grid-points
    ! were distributed according to x-columns. 

    !n_my_grid_points = n_pts_in_block * n_blocks
    !allocate(grid_point_x(n_my_grid_points),grid_point_y(n_my_grid_points),grid_point_z(n_my_grid_points), &
    !     grid_point_block(n_my_grid_points),grid_point_position(n_my_grid_points),STAT=stat)
    !if(stat/=0) call cq_abort("Error allocating grid_point variables: ",n_my_grid_points,stat)
    n_pts = 0
    do nb = 1, n_blocks
       n0_x = (ind_block_x(nb) - 1) * in_block_x
       n0_y = (ind_block_y(nb) - 1) * in_block_y
       n0_z = (ind_block_z(nb) - 1) * in_block_z
       do nz = 1, in_block_z
          do ny = 1, in_block_y
             do nx = 1, in_block_x
                n_pts = n_pts + 1
                grid_point_x(n_pts) = n0_x + nx
                grid_point_y(n_pts) = n0_y + ny
                grid_point_z(n_pts) = n0_z + nz
                grid_point_block(n_pts) = nb
                grid_point_position(n_pts) = nx + (ny - 1) * in_block_x + &
                     (nz - 1) * in_block_y * in_block_x
             enddo
          enddo
       enddo
    enddo
    call stop_timer(tmr_std_indexing)
    return
  end subroutine set_domains
!!***

end module block_module

