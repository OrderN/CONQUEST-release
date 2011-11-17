! $Id$
! -----------------------------------------------------------
! Module blip_grid_transform_module
! -----------------------------------------------------------
! Code area 11: Basis operations
! -----------------------------------------------------------

!!****h* Conquest/blip_grid_transform_module
!!  NAME
!!   blip_grid_transform_module
!!  PURPOSE
!!   Contains all blip to grid transforms: 
!!    blip_to_support_new
!!    blip_to_grad_new
!!    blip_to_gradgrad_new
!!    do_blip_transform_new
!!    do_blip_grad_transform_new
!!    do_blip_gradgrad_transform_new
!!    distribute result
!!   And all inverse transforms:
!!    inverse_blip_transform_new
!!    inverse_blip_to_grad_new
!!    do_inverse_blip_new
!!    do_inverse_blip_to_grad_new
!!    collect_result
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!!   21/6/2000
!!  MODIFICATION HISTORY
!!   17/05/2001 dave
!!    Added ROBODoc header, indented and removed NSF=4 
!!    dependencies
!!   11:29, 12/11/2004 dave 
!!    Removed inappropriate common dependencies
!!   10:14, 2005/07/11 dave
!!    Small changes for intel fortran compiler happiness - moved
!!    declarations of NSF so that they were BEFORE first use of NSF
!!   2006/06/16 07:39 dave
!!    Changing to variable NSF
!!   2008/02/04 08:25 dave
!!    Changing for output to file not stdout
!!   2008/06/10 ast
!!    Added timers
!!   2011/11/15 16:58 dave
!!    Changes throughout to use new blip data structure
!!  SOURCE
module blip_grid_transform_module

 use datatypes
 use global_module, ONLY: iprint_basis, io_lun
 use timer_stdclocks_module, ONLY: start_timer,stop_timer,tmr_std_basis,tmr_std_allocation

 implicit none

  ! RCS tag for object file identification
  character(len=80), save, private :: RCSid = "$Id$"
!!***

contains

!!****f* blip_grid_transform_module/blip_to_support_new *
!!
!!  NAME 
!!   blip_to_support_new 
!!  USAGE
!! 
!!  PURPOSE
!!   Takes blip coefficients and returns values of support 
!!   functions on grid
!!  INPUTS
!! 
!! 
!!  USES
!!   do_blip_transform_new
!!   distribute_result
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!!   June 2000
!!  MODIFICATION HISTORY
!!   17/05/2001 dave
!!    Added ROBODoc header
!!   11/06/2001 dave
!!    Added GenComms
!!   2008/06/10 ast
!!    Added timers
!!   2011/11/15 17:00 dave
!!    Removed blip module and variables passed to transform
!!  SOURCE
!!
  subroutine blip_to_support_new(myid, support)

    use datatypes
    use numbers, ONLY: zero
    use primary_module, ONLY: bundle
    use GenComms, ONLY: my_barrier, cq_abort
    use functions_on_grid, ONLY: gridfunctions
    use species_module, ONLY: nsf_species
    use support_spec_format, ONLY: supports_on_atom
    use timer_module
    use global_module, ONLY: iprint_basis, IPRINT_TIME_THRES2

    implicit none

    ! Passed variables
    integer,intent(in) :: myid, support

    ! Local variables
    integer :: iprim, nsf_send, spec
    integer :: ii
    type(cq_timer) :: tmr_l_tmp1

    call start_timer(tmr_std_basis)
    call start_timer(tmr_l_tmp1,WITH_LEVEL)
    gridfunctions(support)%griddata = zero

    if(bundle%mx_iprim < 1) then
       call cq_abort(' No primary atoms for processor=',myid)
    endif
    do iprim = 1,bundle%mx_iprim
       if(iprint_basis>2.AND.myid==0) write(io_lun,fmt='(2x,"Proc: ",i5," blip-transform for atom ",i5)') myid,iprim
       call my_barrier()
       if( iprim <= bundle%n_prim ) then
          spec = bundle%species(iprim)
          nsf_send = nsf_species(spec)
          call do_blip_transform_new(iprim, supports_on_atom(iprim), nsf_send)
       else
          nsf_send = 0
       endif
       if(iprint_basis>2.AND.myid==0) write(io_lun,fmt='(2x,"Proc: ",i5," distribute result for atom ",i5)') myid,iprim
       call my_barrier()
       call distribute_result(myid,iprim,nsf_send,support)
    end do

    !call time_barrier('after transform',15)

    call stop_print_timer(tmr_l_tmp1,"blip to grid transform",IPRINT_TIME_THRES2)
    call stop_timer(tmr_std_basis)
    return
  end subroutine blip_to_support_new
!!***

!!****f* blip_grid_transform_module/do_blip_transform_new *
!!
!!  NAME 
!!   do_blip_transform_new
!!  USAGE
!! 
!!  PURPOSE
!!   Does the actual transform onto a local grid
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!!   21/6/2000 
!!  MODIFICATION HISTORY
!!   28/6/2000 TM
!!    The rank of data_blip has been changed.
!!   17/05/2001 dave
!!    Added ROBODoc header, indented and removed NSF=4 dependence
!!   2008/06/10 ast
!!    Added timers
!!  SOURCE
!!
  subroutine do_blip_transform_new( iprim, data_blip, nsf)

    use datatypes
    use numbers
    use global_module,   ONLY: rcellx,rcelly,rcellz, area_basis
    use group_module,    ONLY: blocks
    use primary_module,  ONLY: bundle
    use cover_module,    ONLY: BCS_blocks
    use set_blipgrid_module, ONLY: naba_blk_supp
    use comm_array_module,   ONLY: send_array
    use block_module,    ONLY: nx_in_block,ny_in_block,nz_in_block, n_pts_in_block
    use GenComms, ONLY: cq_abort, myid
    use support_spec_format, ONLY: support_function
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use blip, ONLY: blip_info

    ! checking
    !use common, ONLY: inode

    implicit none

    integer,intent(in)     :: iprim
    integer,intent(in)     :: nsf
    type(support_function),intent(in):: data_blip

    ! Local variables
    real(double),allocatable:: inter_1(:,:,:,:)
    real(double),allocatable:: inter_2(:,:,:,:)
    real(double):: rec_sgs, dx, dy, dz, dxxx, dyyy, dzzz
    integer     :: x, y, z, bx, by, bz 
    integer     :: imin,imax, nsf1, spec

    real(double):: ys, yss, dsum(nsf)

    !NEW VARIABLES introduced by TM
    real(double), allocatable, dimension(:) :: splines
    real(double), allocatable, dimension(:,:) :: splines_for_z
    integer     :: ierr
    integer     :: ix_grid,iy_grid,iz_grid,ix,iy,iz
    integer     :: nblkx,nblky,nblkz  ! nos. of integ. pts in a block 
    ! along x,y and z direction
    integer     :: naba_blk,nx_blk,ny_blk,nz_blk,igrid
    real(double) :: dcellx_block,dcelly_block,dcellz_block
    real(double) :: dcellx_grid, dcelly_grid, dcellz_grid

    integer     :: ind_blk,ind_blip
    integer     :: nxmin_grid,nxmax_grid,nymin_grid,nymax_grid,nzmin_grid,nzmax_grid
    integer     :: ncover_yz
    integer     :: stat=0,ii

    integer, allocatable, dimension(:) :: imin_for_z, imax_for_z
    ! inter_1,inter_2,splines,splines_for_z,imin_for_z,imax_for_z
    ! are dynamically allocated variables.
    !**************************************************************************
    !     Start of subroutine

    spec = bundle%species(iprim)
    rec_sgs = one / blip_info(spec)%SupportGridSpacing
    nblkx=nx_in_block
    nblky=ny_in_block
    nblkz=nz_in_block


    dcellx_block=rcellx/blocks%ngcellx; dcellx_grid=dcellx_block/nblkx
    dcelly_block=rcelly/blocks%ngcelly; dcelly_grid=dcelly_block/nblky
    dcellz_block=rcellz/blocks%ngcellz; dcellz_grid=dcellz_block/nblkz

    nxmin_grid=(naba_blk_supp%nxmin(iprim)-1)*nblkx  
    nxmax_grid= naba_blk_supp%nxmax(iprim)*nblkx-1
    nymin_grid=(naba_blk_supp%nymin(iprim)-1)*nblky  
    nymax_grid= naba_blk_supp%nymax(iprim)*nblky-1
    nzmin_grid=(naba_blk_supp%nzmin(iprim)-1)*nblkz  
    nzmax_grid= naba_blk_supp%nzmax(iprim)*nblkz-1
    !check extent
    ierr=0
    if(2*blip_info(spec)%Extent < nxmax_grid-nxmin_grid) then
       write(io_lun,*) myid,' ERROR in do_blip : extent for x dir ::', &
            ' extent nxmin_grid nxmax_grid = ',blip_info(spec)%Extent,nxmin_grid,nxmax_grid
       ierr=1
    endif
    if(2*blip_info(spec)%Extent < nymax_grid-nymin_grid) then
       write(io_lun,*) myid,' ERROR in do_blip : extent for y dir ::', &
            ' extent nymin_grid nymax_grid = ',blip_info(spec)%Extent,nymin_grid,nymax_grid
       ierr=1
    endif
    if(2*blip_info(spec)%Extent < nzmax_grid-nzmin_grid) then
       write(io_lun,*) myid,' ERROR in do_blip : extent for z dir ::', &
            ' extent nymin_grid nymax_grid = ',blip_info(spec)%Extent,nzmin_grid,nzmax_grid
       ierr=1
    endif
    if(ierr /= 0) call cq_abort('do_blip ',ierr)

    call start_timer(tmr_std_allocation)
    allocate(inter_1(nsf,-blip_info(spec)%BlipArraySize:blip_info(spec)%BlipArraySize, &
         -blip_info(spec)%BlipArraySize:blip_info(spec)%BlipArraySize,2*blip_info(spec)%Extent+1),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to inter_1 !',stat)
    call reg_alloc_mem(area_basis,nsf*(2*blip_info(spec)%BlipArraySize+1)*(2*blip_info(spec)%BlipArraySize+1) &
         *(2*blip_info(spec)%Extent+1),type_dbl)
    call stop_timer(tmr_std_allocation)
    inter_1=zero
    call start_timer(tmr_std_allocation)
    allocate(inter_2(nsf,-blip_info(spec)%BlipArraySize:blip_info(spec)%BlipArraySize, &
         2*blip_info(spec)%Extent+1 ,2*blip_info(spec)%Extent+1),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to inter_2 !',blip_info(spec)%BlipArraySize,blip_info(spec)%Extent)
    call reg_alloc_mem(area_basis,nsf*(2*blip_info(spec)%BlipArraySize+1)*&
         (2*blip_info(spec)%Extent+1)*(2*blip_info(spec)%Extent+1),type_dbl)
    allocate(splines( -blip_info(spec)%BlipArraySize:blip_info(spec)%BlipArraySize ), &
         splines_for_z(-blip_info(spec)%BlipArraySize:blip_info(spec)%BlipArraySize,1:2*blip_info(spec)%Extent+1))
    call reg_alloc_mem(area_basis,(2*blip_info(spec)%BlipArraySize+1)*(2*blip_info(spec)%Extent+2),type_dbl)
    allocate(imin_for_z(2*blip_info(spec)%Extent+1),imax_for_z(2*blip_info(spec)%Extent+1),STAT=stat)
    call reg_alloc_mem(area_basis,2*(2*blip_info(spec)%Extent+1),type_dbl)
    call stop_timer(tmr_std_allocation)
    inter_2=zero
    imin_for_z = 0
    imax_for_z = 0
    splines = zero
    splines_for_z = zero

    !x,y,z are used for integration grids
    !bx,by,bz are used for support grids
    DO x = nxmin_grid, nxmax_grid
       ix_grid=x-nxmin_grid+1
       dxxx=x*dcellx_grid-bundle%xprim(iprim) ! = dxx+dx in old version
       imin = MAX(-blip_info(spec)%BlipArraySize,int((dxxx-blip_info(spec)%BlipWidth*half)*rec_sgs))
       imax = MIN( blip_info(spec)%BlipArraySize,int((dxxx+blip_info(spec)%BlipWidth*half)*rec_sgs))

       if(imin >  blip_info(spec)%BlipArraySize) imin= blip_info(spec)%BlipArraySize
       if(imax < -blip_info(spec)%BlipArraySize) imax=-blip_info(spec)%BlipArraySize

       !OLD dx = rionxi*n_grid_x
       !OLD dx = (- dx + anint(dx))*x_grid*r_super_x
       !OLD dxx = x * r_super_x * x_grid
       !NEW dx+dxx=dxxx
       ! I (TM) think the following expresseion is enough for imin, 
       ! but the difference of the efficiency is small.
       !  imin = MAX(-blip_info(spec)%BlipArraySize,
       !              int((dxx+dx-blip_info(spec)%BlipWidth*half)*rec_sgs+one-very_small))
       DO bx = imin, imax
          yss = dxxx - bx * blip_info(spec)%SupportGridSpacing 
          ys = ABS( yss )
          ys = ys * blip_info(spec)%FourOnBlipWidth
          IF( ys .LE. one ) THEN
             splines( bx ) = one + &
                  ys * ys * ( 0.75_double * ys - 1.5_double )
          ELSE IF( ys .LE. two ) THEN
             splines( bx ) = 0.25_double * ( two - ys ) *  &
                  ( two - ys ) * ( two - ys ) 
          ELSE
             splines( bx ) = zero
          END IF
       END DO
       DO bz = -blip_info(spec)%BlipArraySize, blip_info(spec)%BlipArraySize
          DO by = -blip_info(spec)%BlipArraySize, blip_info(spec)%BlipArraySize
             dsum = zero
             DO bx = imin,imax
                ind_blip = blip_info(spec)%blip_number(bx,by,bz)
                IF( ind_blip .NE. 0 ) THEN
                   ys = splines( bx )
                   do nsf1 = 1,nsf
                      dsum(nsf1) = dsum(nsf1) + ys * data_blip%supp_func(nsf1)%coefficients(ind_blip)
                   enddo
                END IF
             END DO
             do nsf1 = 1,nsf
                inter_1( nsf1,by,bz, ix_grid ) = dsum(nsf1)
             enddo
          END DO
       END DO
    END DO

    DO y = nymin_grid, nymax_grid
       iy_grid=y-nymin_grid+1
       dyyy=y*dcelly_grid-bundle%yprim(iprim)
       imin = MAX(-blip_info(spec)%BlipArraySize,int((dyyy-blip_info(spec)%BlipWidth*half)*rec_sgs))
       imax = MIN( blip_info(spec)%BlipArraySize,int((dyyy+blip_info(spec)%BlipWidth*half)*rec_sgs))
       if(imin >  blip_info(spec)%BlipArraySize) imin= blip_info(spec)%BlipArraySize
       if(imax < -blip_info(spec)%BlipArraySize) imax=-blip_info(spec)%BlipArraySize

       DO by = imin,imax
          yss = dyyy - by * blip_info(spec)%SupportGridSpacing 
          ys = ABS( yss )
          ys = ys * blip_info(spec)%FourOnBlipWidth
          IF( ys .LE. one ) THEN
             splines( by ) = one +  &
                  ys * ys * ( 0.75_double * ys - 1.5_double )
          ELSE IF( ys .LE. two ) THEN
             splines( by ) = 0.25_double * ( two - ys ) *  &
                  ( two - ys ) * ( two - ys )  
          ELSE
             splines( by ) = zero
          END IF
       END DO
       DO x = nxmin_grid, nxmax_grid
          ix_grid=x-nxmin_grid+1
          DO bz = -blip_info(spec)%BlipArraySize, blip_info(spec)%BlipArraySize
             dsum = zero
             DO by = imin,imax
                !DO by = -blip_info(spec)%BlipArraySize,blip_info(spec)%BlipArraySize
                ys = splines( by )
                do nsf1 = 1,nsf
                   dsum(nsf1) = dsum(nsf1)+ys*inter_1( nsf1, by, bz, ix_grid)
                enddo
             END DO
             do nsf1 = 1,nsf
                inter_2( nsf1, bz,ix_grid,iy_grid) = dsum(nsf1)
             enddo
          END DO
       END DO
    END DO

    splines_for_z = BIG

    DO z = nzmin_grid, nzmax_grid
       iz_grid=z-nzmin_grid+1
       !dzz = z * r_super_z * z_grid
       dzzz=z*dcellz_grid-bundle%zprim(iprim)
       imin_for_z(iz_grid) = MAX(-blip_info(spec)%BlipArraySize,&
            int((dzzz-blip_info(spec)%BlipWidth*half)*rec_sgs))
       imax_for_z(iz_grid) = MIN( blip_info(spec)%BlipArraySize,&
            int((dzzz+blip_info(spec)%BlipWidth*half)*rec_sgs))

       if(imin_for_z(iz_grid) >  blip_info(spec)%BlipArraySize) then
          imin_for_z(iz_grid)= blip_info(spec)%BlipArraySize
       endif
       if(imax_for_z(iz_grid) < -blip_info(spec)%BlipArraySize) then
          imax_for_z(iz_grid)=-blip_info(spec)%BlipArraySize
       endif

       DO bz = imin_for_z(iz_grid),imax_for_z(iz_grid)
          yss = dzzz - bz * blip_info(spec)%SupportGridSpacing 
          ys = ABS( yss )
          ys = ys * blip_info(spec)%FourOnBlipWidth
          IF( ys .LE. one ) THEN
             splines_for_z(bz,iz_grid) = one +                      &
                  ys * ys * ( 0.75_double * ys - 1.5_double )
          ELSE IF( ys .LE. two ) THEN
             splines_for_z(bz,iz_grid) = 0.25_double * ( two - ys ) *    &
                  ( two - ys ) * ( two - ys )  
          ELSE
             splines_for_z(bz,iz_grid) = zero
          END IF
       END DO
    ENDDO

    ncover_yz=BCS_blocks%ncovery*BCS_blocks%ncoverz
    igrid=0

    call start_timer(tmr_std_allocation)
    allocate(send_array(naba_blk_supp%no_naba_blk(iprim)*nsf*n_pts_in_block),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating send_array in do_blip_transform: ",&
         naba_blk_supp%no_naba_blk(iprim)*nsf*n_pts_in_block)
    call reg_alloc_mem(area_basis,naba_blk_supp%no_naba_blk(iprim)*nsf*n_pts_in_block,type_dbl)
    call stop_timer(tmr_std_allocation)
    send_array(:) = zero
    DO naba_blk=1,naba_blk_supp%no_naba_blk(iprim)! naba blks in NOPG order
       ! iprim : primary seq. no. of the atom
       ind_blk=naba_blk_supp%list_naba_blk(naba_blk,iprim) !CC in BCS_blocks
       ind_blk=ind_blk-1

       nx_blk=ind_blk/ncover_yz
       ny_blk=(ind_blk-nx_blk*ncover_yz)/BCS_blocks%ncoverz
       nz_blk=ind_blk-nx_blk*ncover_yz-ny_blk*BCS_blocks%ncoverz

       nx_blk=nx_blk-BCS_blocks%nspanlx+BCS_blocks%nx_origin
       ny_blk=ny_blk-BCS_blocks%nspanly+BCS_blocks%ny_origin
       nz_blk=nz_blk-BCS_blocks%nspanlz+BCS_blocks%nz_origin

       do iz=1,nblkz    ! z-direction in a block
          z=nblkz*(nz_blk-1)+iz-1
          iz_grid=z-nzmin_grid+1
          if(iz_grid > 2*blip_info(spec)%Extent+1) then
             call cq_abort(' ERROR in do_blip_transform_new iz_grid = ', iz_grid,2*blip_info(spec)%Extent+1)
          endif

          imin=imin_for_z(iz_grid)
          imax=imax_for_z(iz_grid)

          do iy=1,nblky    ! y-direction in a block
             y=nblky*(ny_blk-1)+iy-1
             iy_grid=y-nymin_grid+1
             if(iy_grid > 2*blip_info(spec)%Extent+1) then
                call cq_abort(' ERROR in do_blip_transform_new iy_grid = ', iy_grid,2*blip_info(spec)%Extent+1)
             endif

             do ix=1,nblkx    ! x-direction in a block
                x=nblkx*(nx_blk-1)+ix-1
                ix_grid=x-nxmin_grid+1
                if(ix_grid > 2*blip_info(spec)%Extent+1) then
                   call cq_abort(' ERROR in do_blip_transform_new ix_grid = ', ix_grid,2*blip_info(spec)%Extent+1)
                endif

                igrid=igrid+1  ! seq. no. of integration grids
                dsum = zero
                DO bz = imin,imax
                   !DO bz = -blip_info(spec)%BlipArraySize,blip_info(spec)%BlipArraySize
                   ys = splines_for_z(bz,iz_grid)
                   do nsf1 = 1,nsf
                      dsum(nsf1)=dsum(nsf1)+ys*inter_2(nsf1,bz,ix_grid,iy_grid)
                   enddo
                END DO
                do nsf1 = 1,nsf
                   send_array(nsf*(igrid-1)+nsf1) = dsum(nsf1)
                enddo
             enddo    ! x-direction in a block
          enddo    ! y-direction in a block
       enddo    ! z-direction in a block
    END DO    ! Loop over naba blocks for this primary atom
    call start_timer(tmr_std_allocation)
    deallocate(inter_1,inter_2,STAT=stat)
    if(stat/=0) call cq_abort('Error deallocating memory to inter_1 !',stat)
    call reg_dealloc_mem(area_basis,nsf*(2*blip_info(spec)%BlipArraySize+1)*(2*blip_info(spec)%BlipArraySize+1)*&
         (2*blip_info(spec)%Extent+1),type_dbl)
    call reg_dealloc_mem(area_basis,nsf*(2*blip_info(spec)%BlipArraySize+1)*&
         (2*blip_info(spec)%Extent+1)*(2*blip_info(spec)%Extent+1),type_dbl)
    deallocate(splines,splines_for_z,STAT=stat)
    call reg_dealloc_mem(area_basis,(2*blip_info(spec)%BlipArraySize+1)*(2*blip_info(spec)%Extent+2),type_dbl)
    if(stat/=0) call cq_abort('Error deallocating memory to splines !',stat)
    deallocate(imin_for_z,imax_for_z,STAT=stat)
    call reg_dealloc_mem(area_basis,2*(2*blip_info(spec)%Extent+1),type_dbl)
    call stop_timer(tmr_std_allocation)
    return
  end subroutine do_blip_transform_new
!!***

 !---------------------------------------------------------------
 !sbrt distribute_result
 !  this subroutine distributes the blip-grid transformed 
 !  values on grid points to their domain responsible node.
 !       
 !---------------------------------------------------------------
  subroutine distribute_result(myid,iprim,nsf_send,support)
    !  new version 
    !     21/6/2000 Tsuyoshi Miyazaki
    !  if we include nsf to dummy arguments, we can use this subroutine
    !  for the case where the number of support functions depends on
    !  each atom. (Of course, it is needed to make send_array correctly
    !  in do_blip_transform.)
    !   26/6/2000 Now we assume that this subroutine is called once to 
    !             send blip-grid transformed support functions of one
    !             primary atom. If the number of support func. is
    !             different for different atoms, nunit and nsize varies
    !             at each call.
    !
    use datatypes
    use numbers, ONLY : zero
    use mpi
    use set_blipgrid_module, ONLY: comBG,naba_atm
    use comm_array_module,   ONLY: send_array,recv_array
    use block_module,        ONLY: n_pts_in_block  ! = blocks%nm_group(:)
    use GenComms, ONLY: my_barrier, cq_abort
    use functions_on_grid, ONLY: gridfunctions, fn_on_grid
    use species_module, ONLY: nsf_species
    use global_module, ONLY: species_glob, sf, area_basis
    use atoms, ONLY: atoms_on_node
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none
    integer,intent(in) :: myid,iprim,nsf_send
    integer :: support

    integer :: mynode, nunit_send, nunit_recv, ind_alp_i_blk
    integer :: tag
    integer :: nsend_req(comBG%mx_send_node),ierr_send(comBG%mx_send_node)
    integer :: nrecv_stat(MPI_STATUS_SIZE,comBG%mx_recv_node),&
         ierr_recv(comBG%mx_recv_node)
    !TM 20/12/2000
    integer :: nwait_stat(MPI_STATUS_SIZE),ierr,irc
    integer :: inode,jnode,nnd_rem,ibegin,nsize,my_ibegin, nsf_recv
    integer :: igrid,isf,prim_blk,naba_atm_tmp,loc1,loc2,ind1,ind2, stat
    !integer,save :: ipair
    integer :: ipair
    integer :: isend
    real(double),pointer :: recv_ptr(:)

    mynode=myid+1
    ! nsf is set by the calling routine
    tag=1
    !send blip-grid transformed support functions
    !  if recv_node is my node, it just keeps the first address of the
    !  corresponding data in send_array.
    !  

    if(comBG%no_recv_node(iprim) > 0) then
       nunit_send=nsf_send*n_pts_in_block   
       do inode=1,comBG%no_recv_node(iprim)
          nnd_rem=comBG%list_recv_node(inode,iprim)
          if(inode == 1) then
             ibegin=1
          else
             ibegin=ibegin+comBG%no_naba_blk(inode-1,iprim)*nunit_send
          endif

          nsize=comBG%no_naba_blk(inode,iprim)*nunit_send
          if(nnd_rem == mynode) then
             my_ibegin=ibegin
          else
             call MPI_issend(send_array(ibegin),nsize,MPI_DOUBLE_PRECISION,nnd_rem-1, &
                  tag,MPI_COMM_WORLD,nsend_req(inode),ierr_send(inode))
             if(ierr_send(inode) /= 0) then
                call cq_abort('ERROR in MPI_issend for ',ierr_send(inode),inode)
             endif
          endif
       enddo
    endif ! if there are receiving nodes

    call my_barrier()
    !receive

    if(comBG%no_send_node(iprim) > 0) then
       isend=comBG%ibeg_recv_call(iprim)-1
       do jnode=1,comBG%no_send_node(iprim)
          isend=isend+1
          nnd_rem=comBG%list_send_node(jnode,iprim)
          nsf_recv = nsf_species(species_glob(atoms_on_node(iprim,nnd_rem)))
          nunit_recv = nsf_recv*n_pts_in_block
          nsize=comBG%no_sent_pairs(jnode,iprim)*nunit_recv
          call start_timer(tmr_std_allocation)
          allocate(recv_array(nsize),STAT=stat)
          if(stat/=0) call cq_abort("Error allocating recv_array in distribute_result: ",nsize)
          call reg_alloc_mem(area_basis,nsize,type_dbl)
          call stop_timer(tmr_std_allocation)
          if(nnd_rem == mynode) then
             if(nsize < 1) call cq_abort('ERROR recv_size in distribute_result',nsize)

             recv_ptr => send_array(my_ibegin:my_ibegin+nsize-1)
          else
             if(nsize < 1) call cq_abort('ERROR recv_size in distribute_result',nsize)

             call MPI_recv(recv_array,nsize,MPI_DOUBLE_PRECISION,nnd_rem-1, &
                  tag,MPI_COMM_WORLD,nrecv_stat(:,jnode),ierr_recv(jnode))
             if(ierr_recv(jnode) /= 0) then
                call cq_abort('ERROR in MPI_recv',ierr_recv(jnode),jnode)
             endif
             recv_ptr => recv_array(1:nsize)
          endif

          do ipair=1,comBG%no_sent_pairs(jnode,iprim) !naba blks of the iprim-th atm
             prim_blk=comBG%table_blk(ipair,isend)     ! primary block in my domain
             naba_atm_tmp=comBG%table_pair(ipair,isend)! naba atm of the primary blk
             loc1= nunit_recv*(ipair-1)
             ind_alp_i_blk = naba_atm(sf)%ibegin_blk_orb(prim_blk)+ &
                  naba_atm(sf)%ibeg_orb_atom(naba_atm_tmp, prim_blk) -1 
             ! naba_atm(sf)%ibegin_blk_orb(iprim_blk): initial position of support for the present primary block.
             ! naba_atm(sf)%ibeg_orb_atom(naba, iprim_blk)
             ! : initial position of the present naba atom (naba) in the neighbour-atom orbitals of the 
             ! present primary block.
             !loc2 = n_pts_in_block * (ind_alp_i_blk-1) + 1
             loc2 = n_pts_in_block * (ind_alp_i_blk-1)

             if(loc1+nunit_recv>nsize) call cq_abort('ERROR in distribute_result, loc1: ',loc1+nunit_recv,nsize)
             if(loc2+nunit_recv>gridfunctions(support)%size) &
                  call cq_abort('ERROR in distribute_result :loc2',loc2+nunit_recv,gridfunctions(support)%size)

             !!------ TRANSPOSE(grid<->orbital(nsf) ) within this pair -----
             !!   We can NOT do this in making send_array as is done in 
             !!   the old version. In the old version, we use two arrays to
             !!   do blip-grid transforms and distribute the results.
             !!   The size of the former array is about twice (the volume of 
             !!   cube/sphere) of the one of present send_array.
             !!   Now send_array is arranged like (nsf,naba_blks) and it can 
             !!   be used directly in communications because naba_blks 
             !!   are arranged in NOPG(node order periodically grouped) order.
             !!
             !!   We can change the order of the following two do-loops.
             !    ind1=loc1
             ! do igrid=1,n_pts_in_block
             !  do isf=1,nsf_send
             !    ind1=ind1+1  !  <=>  ind1=(igrid-1)*nsf_send+isf + loc1
             !    ind2=(isf-1)*n_pts_in_block+igrid
             !    support(loc2+ind2)=recv_ptr(ind1)
             !  enddo
             ! enddo
             do isf=1,nsf_recv
                ind2=(isf-1)*n_pts_in_block+loc2
                do igrid=1,n_pts_in_block
                   ind1=(igrid-1)*nsf_recv+isf 
                   ind2=ind2+1
                   gridfunctions(support)%griddata(ind2)=recv_ptr(ind1+loc1)
                enddo
             enddo
             !!------ TRANSPOSE

          enddo ! loop over pairs of block and atoms
          call start_timer(tmr_std_allocation)
          deallocate(recv_array,STAT=stat)
          if(stat/=0) call cq_abort("Error deallocating recv_array in distribute_result: ",nsize)
          call reg_dealloc_mem(area_basis,nsize,type_dbl)
          call stop_timer(tmr_std_allocation)
       enddo ! loop over sending nodes
    endif ! if there are sending nodes

    !Check whether we have finished all MPI_issend
    !     20/12/2000 Tsuyoshi Miyazaki 
    if(comBG%no_recv_node(iprim) > 0) then
       do inode=1,comBG%no_recv_node(iprim)
          nnd_rem=comBG%list_recv_node(inode,iprim)
          if(nnd_rem /= mynode) then
             call MPI_WAIT(nsend_req(inode),nwait_stat,ierr)
             if(ierr /= 0) call cq_abort('ERROR in calling MPI_WAIT in distribute_result',ierr)
          endif      ! if it is remote...
       enddo      ! Loop over remote nodes
    endif       ! if there are sending nodes for iprim...

    call my_barrier()
    if(allocated(send_array)) then
       call start_timer(tmr_std_allocation)
       call reg_dealloc_mem(area_basis,size(send_array),type_dbl)
       deallocate(send_array,STAT=stat)
       if(stat/=0) call cq_abort("Error deallocating send_array in do_blip_transform: ",stat)
       call stop_timer(tmr_std_allocation)
    end if
    return
  end subroutine distribute_result

!!****f* blip_grid_transform_module/blip_to_grad_new *
!!
!!  NAME 
!!   blip_to_grad_new
!!  USAGE
!! 
!!  PURPOSE
!!   Takes the blip coefficients and returns their derivative on
!!   the grid
!!  INPUTS
!! 
!! 
!!  USES
!!   do_blip_grad_transform_new
!!   distribute_result
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!!   28/12/2000
!!  MODIFICATION HISTORY
!!   17/05/2001 dave
!!    Added ROBODoc header, removed unnecessary passed arguments
!!   11/06/2001 dave
!!    Added GenComms
!!   2008/06/10 ast
!!    Added timers
!!   2011/11/15 17:00 dave
!!    Removed blip module and variables passed to transform
!!  SOURCE
!!
  subroutine blip_to_grad_new(myid, direction, support)

    use datatypes
    use numbers
    use blip, ONLY: blip_info
    use primary_module, ONLY: bundle
    use GenComms, ONLY: my_barrier, cq_abort
    use functions_on_grid, ONLY: gridfunctions
    use species_module, ONLY: nsf_species
    use support_spec_format, ONLY: supports_on_atom

    implicit none

    ! Passed variables
    integer,intent(in) :: myid, direction, support

    ! Local variables
    integer :: iprim, nsf_send, spec
    integer :: ii

    call start_timer(tmr_std_basis)
    gridfunctions(support)%griddata = zero

    if(bundle%mx_iprim < 1) call cq_abort('ERROR in blip_to_grad_new: no primary atoms')
    do iprim = 1,bundle%mx_iprim
       call my_barrier()
       if( iprim <= bundle%n_prim ) then
          spec = bundle%species(iprim)
          nsf_send = nsf_species(spec)
          call do_blip_grad_transform_new(direction, iprim, supports_on_atom(iprim), nsf_send)
       else
          nsf_send = 0
       endif
       call distribute_result(myid,iprim,nsf_send,support)
    end do
    call stop_timer(tmr_std_basis)
    return
  end subroutine blip_to_grad_new
!!***

!!****f* blip_grid_transform_module/do_blip_grad_transform_new *
!!
!!  NAME 
!!   do_blip_grad_transform_new
!!  USAGE
!! 
!!  PURPOSE
!!   Does the blip to grad transform onto a local grid
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!!   28/12/2000
!!  MODIFICATION HISTORY
!!   17/05/2001 dave
!!    Added ROBODoc header, indented and removed NSF=4 dependency
!!   2008/06/10 ast
!!    Added timers
!!  SOURCE
!!
  subroutine do_blip_grad_transform_new(direction, iprim, data_blip, nsf)

    use datatypes
    use numbers
    use global_module,   ONLY: rcellx,rcelly,rcellz, area_basis
    use group_module,    ONLY: blocks
    use primary_module,  ONLY: bundle
    use cover_module,    ONLY: BCS_blocks
    use set_blipgrid_module, ONLY: naba_blk_supp
    use comm_array_module,   ONLY: send_array
    use block_module,    ONLY: nx_in_block,ny_in_block,nz_in_block, n_pts_in_block
    use GenComms,        ONLY: cq_abort
    use support_spec_format, ONLY: support_function
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use blip, ONLY: blip_info

    implicit none

    integer,intent(in)     :: direction
    integer,intent(in)     :: iprim
    integer,intent(in)     :: NSF
    type(support_function),intent(in):: data_blip

    ! Local variables
    real(double),allocatable:: inter_1(:,:,:,:)
    real(double),allocatable:: inter_2(:,:,:,:)
    real(double):: rec_sgs, dx, dy, dz, dxxx, dyyy, dzzz
    integer     :: x, y, z, bx, by, bz 
    integer     :: imin,imax, nsf1, spec

    real(double):: ys, yss, dsum(NSF)

    !NEW VARIABLES introduced by TM
    real(double), allocatable, dimension(:) :: splines
    real(double), allocatable, dimension(:,:) :: splines_for_z
    integer, allocatable, dimension(:) :: imin_for_z, imax_for_z
    integer     :: ierr
    integer     :: ix_grid,iy_grid,iz_grid,ix,iy,iz
    integer     :: nblkx,nblky,nblkz  ! nos. of integ. pts in a block 
    ! along x,y and z direction
    integer     :: naba_blk,nx_blk,ny_blk,nz_blk,igrid
    real(double) :: dcellx_block,dcelly_block,dcellz_block
    real(double) :: dcellx_grid, dcelly_grid, dcellz_grid

    integer     :: ind_blk,ind_blip
    integer     :: nxmin_grid,nxmax_grid,nymin_grid,nymax_grid,&
         nzmin_grid,nzmax_grid
    integer     :: ncover_yz
    integer     :: stat=0,ii

    !real(double) :: sum_send_array1,sum_send_array2,&
         !                sum_send_array3,sum_send_array4
    ! inter_1,inter_2,splines,splines_for_z,imin_for_z,imax_for_z
    ! are dynamically allocated variables.

    spec = bundle%species(iprim)
    rec_sgs = one / blip_info(spec)%SupportGridSpacing
    nblkx=nx_in_block
    nblky=ny_in_block
    nblkz=nz_in_block

    dcellx_block=rcellx/blocks%ngcellx; dcellx_grid=dcellx_block/nblkx
    dcelly_block=rcelly/blocks%ngcelly; dcelly_grid=dcelly_block/nblky
    dcellz_block=rcellz/blocks%ngcellz; dcellz_grid=dcellz_block/nblkz

    nxmin_grid=(naba_blk_supp%nxmin(iprim)-1)*nblkx  
    nxmax_grid= naba_blk_supp%nxmax(iprim)*nblkx-1
    nymin_grid=(naba_blk_supp%nymin(iprim)-1)*nblky  
    nymax_grid= naba_blk_supp%nymax(iprim)*nblky-1
    nzmin_grid=(naba_blk_supp%nzmin(iprim)-1)*nblkz  
    nzmax_grid= naba_blk_supp%nzmax(iprim)*nblkz-1
    !check extent & MAXBAS.
    ierr=0
    if(2*blip_info(spec)%Extent < nxmax_grid-nxmin_grid) then
       write(io_lun,*) ' ERROR in do_blip_transform : extent for x direction ::', &
            ' extent nxmin_grid nxmax_grid = ',blip_info(spec)%Extent,nxmin_grid,nxmax_grid
       ierr=1
    endif
    if(2*blip_info(spec)%Extent < nymax_grid-nymin_grid) then
       write(io_lun,*) ' ERROR in do_blip_transform : extent for y direction ::', &
            ' extent nymin_grid nymax_grid = ',blip_info(spec)%Extent,nymin_grid,nymax_grid
       ierr=1
    endif
    if(2*blip_info(spec)%Extent < nzmax_grid-nzmin_grid) then
       write(io_lun,*) ' ERROR in do_blip_transform : extent for z direction ::', &
            ' extent nymin_grid nymax_grid = ',blip_info(spec)%Extent,nzmin_grid,nzmax_grid
       ierr=1
    endif
    if(ierr /= 0) call cq_abort('ERROR in do_blip_grad_transform')
    call start_timer(tmr_std_allocation)
    allocate(inter_1(NSF,-blip_info(spec)%BlipArraySize:blip_info(spec)%BlipArraySize,&
         -blip_info(spec)%BlipArraySize:blip_info(spec)%BlipArraySize,2*blip_info(spec)%Extent+1),STAT=stat)
    if(stat/=0) call cq_abort('ERROR in do_blip_grad_transform: allocation of inter_1')
    call reg_alloc_mem(area_basis,nsf*(2*blip_info(spec)%BlipArraySize+1)*(2*blip_info(spec)%BlipArraySize+1)*&
         (2*blip_info(spec)%Extent+1),type_dbl)
    call stop_timer(tmr_std_allocation)
    inter_1=zero
    call start_timer(tmr_std_allocation)
    allocate(inter_2(NSF,-blip_info(spec)%BlipArraySize:blip_info(spec)%BlipArraySize,&
         2*blip_info(spec)%Extent+1,2*blip_info(spec)%Extent+1),STAT=stat)
    if(stat/=0) call cq_abort('ERROR in do_blip_grad_transform: allocation of inter_2')
    call reg_alloc_mem(area_basis,nsf*(2*blip_info(spec)%BlipArraySize+1)*&
         (2*blip_info(spec)%Extent+1)*(2*blip_info(spec)%Extent+1),type_dbl)
    allocate(splines( -blip_info(spec)%BlipArraySize:blip_info(spec)%BlipArraySize ), &
         splines_for_z(-blip_info(spec)%BlipArraySize:blip_info(spec)%BlipArraySize,1:2*blip_info(spec)%Extent+1))
    call reg_alloc_mem(area_basis,(2*blip_info(spec)%BlipArraySize+1)*(2*blip_info(spec)%Extent+2),type_dbl)
    allocate(imin_for_z(2*blip_info(spec)%Extent+1),imax_for_z(2*blip_info(spec)%Extent+1),STAT=stat)
    call reg_alloc_mem(area_basis,2*(2*blip_info(spec)%Extent+1),type_dbl)
    call stop_timer(tmr_std_allocation)
    inter_2=zero
    imin_for_z = 0
    imax_for_z = 0
    splines = zero
    splines_for_z = zero

    !x,y,z are used for integration grids
    !bx,by,bz are used for support grids
    DO x = nxmin_grid, nxmax_grid
       ix_grid=x-nxmin_grid+1
       dxxx=x*dcellx_grid-bundle%xprim(iprim) ! = dxx+dx in old version
       imin = MAX(-blip_info(spec)%BlipArraySize,int((dxxx-blip_info(spec)%BlipWidth*half)*rec_sgs))
       imax = MIN( blip_info(spec)%BlipArraySize,int((dxxx+blip_info(spec)%BlipWidth*half)*rec_sgs))

       !OLD dx = rionxi*n_grid_x
       !OLD dx = (- dx + anint(dx))*x_grid*r_super_x
       !OLD dxx = x * r_super_x * x_grid
       !NEW dx+dxx=dxxx
       ! I (TM) think the following expresseion is enough for imin, 
       ! but the difference of the efficiency is small.
       !  imin = MAX(-blip_info(spec)%BlipArraySize,
       !              int((dxx+dx-blip_info(spec)%BlipWidth*half)*rec_sgs+one-very_small))
       DO bx = imin, imax
          yss = dxxx - bx * blip_info(spec)%SupportGridSpacing 
          ys = ABS( yss )
          ys = ys * blip_info(spec)%FourOnBlipWidth
          if (yss.ge.zero) then
             yss = blip_info(spec)%FourOnBlipWidth
          else
             yss = -blip_info(spec)%FourOnBlipWidth
          end if
          IF( ys .LE. one ) THEN
             if (direction.ne.1) then
                splines( bx ) = one + &
                     ys * ys * ( 0.75_double * ys - 1.5_double )
             else
                splines( bx ) = yss * ys * (2.25_double * ys - 3.0_double )
             end if

          ELSE IF( ys .LE. two ) THEN
             if (direction.ne.1) then
                splines( bx ) = 0.25_double * ( two - ys ) *  &
                     ( two - ys ) * ( two - ys ) 
             else
                splines( bx ) = -yss * 0.75_double * ( two - ys ) * ( two - ys )
             end if

          ELSE
             splines( bx ) = zero
          END IF
       END DO
       DO bz = -blip_info(spec)%BlipArraySize, blip_info(spec)%BlipArraySize
          DO by = -blip_info(spec)%BlipArraySize, blip_info(spec)%BlipArraySize
             dsum = zero
             DO bx = imin,imax
                ind_blip = blip_info(spec)%blip_number(bx,by,bz)
                IF( ind_blip .NE. 0 ) THEN
                   ys = splines( bx )
                   do nsf1 = 1,NSF
                      dsum(nsf1)=dsum(nsf1)+ys*data_blip%supp_func(nsf1)%coefficients(ind_blip)
                   enddo
                END IF
             END DO
             do nsf1 = 1,NSF
                inter_1( nsf1,by,bz, ix_grid ) = dsum(nsf1)
             enddo
          END DO
       END DO
    END DO

    DO y = nymin_grid, nymax_grid
       iy_grid=y-nymin_grid+1
       !dyy = y * r_super_y * y_grid
       dyyy=y*dcelly_grid-bundle%yprim(iprim)
       imin = MAX(-blip_info(spec)%BlipArraySize,int((dyyy-blip_info(spec)%BlipWidth*half)*rec_sgs))
       imax = MIN( blip_info(spec)%BlipArraySize,int((dyyy+blip_info(spec)%BlipWidth*half)*rec_sgs))
       DO by = imin,imax
          yss = dyyy - by * blip_info(spec)%SupportGridSpacing 
          ys = ABS( yss )
          ys = ys * blip_info(spec)%FourOnBlipWidth
          if (yss.ge.zero) then
             yss = blip_info(spec)%FourOnBlipWidth
          else
             yss = -blip_info(spec)%FourOnBlipWidth
          end if
          IF( ys .LE. one ) THEN
             if (direction.ne.2) then
                splines( by ) = one +  &
                     ys * ys * ( 0.75_double * ys - 1.5_double )
             else
                splines( by ) = yss * ys * (2.25_double * ys - 3.0_double )
             end if
          ELSE IF( ys .LE. two ) THEN
             if (direction.ne.2) then
                splines( by ) = 0.25_double * ( two - ys ) *  &
                     ( two - ys ) * ( two - ys )  
             else
                splines( by ) = -yss * 0.75_double * ( two - ys ) * ( two - ys )
             end if
          ELSE
             splines( by ) = zero
          END IF
       END DO
       DO x = nxmin_grid, nxmax_grid
          ix_grid=x-nxmin_grid+1
          DO bz = -blip_info(spec)%BlipArraySize, blip_info(spec)%BlipArraySize
             dsum = zero
             DO by = imin,imax
                ys = splines( by )
                do nsf1 = 1,NSF
                   dsum(nsf1)=dsum(nsf1)+ys*inter_1( nsf1, by, bz, ix_grid)
                enddo
             END DO
             do nsf1 = 1,NSF
                inter_2( nsf1, bz,ix_grid,iy_grid) = dsum(nsf1)
             enddo
          END DO
       END DO
    END DO

    DO z = nzmin_grid, nzmax_grid
       iz_grid=z-nzmin_grid+1
       !dzz = z * r_super_z * z_grid
       dzzz=z*dcellz_grid-bundle%zprim(iprim)
       imin_for_z(iz_grid) = MAX(-blip_info(spec)%BlipArraySize,&
            int((dzzz-blip_info(spec)%BlipWidth*half)*rec_sgs))
       imax_for_z(iz_grid) = MIN( blip_info(spec)%BlipArraySize,&
            int((dzzz+blip_info(spec)%BlipWidth*half)*rec_sgs))
       DO bz = imin_for_z(iz_grid),imax_for_z(iz_grid)
          yss = dzzz - bz * blip_info(spec)%SupportGridSpacing 
          ys = ABS( yss )
          ys = ys * blip_info(spec)%FourOnBlipWidth
          if (yss.ge.zero) then
             yss = blip_info(spec)%FourOnBlipWidth
          else
             yss = -blip_info(spec)%FourOnBlipWidth
          end if
          IF( ys .LE. one ) THEN
             if (direction.ne.3) then
                splines_for_z(bz,iz_grid) = one +                      &
                     ys * ys * ( 0.75_double * ys - 1.5_double )
             else
                splines_for_z(bz,iz_grid) = yss * ys * (2.25_double * ys - 3.0_double )
             end if

          ELSE IF( ys .LE. two ) THEN
             if (direction.ne.3) then
                splines_for_z(bz,iz_grid) = 0.25_double * ( two - ys ) *    &
                     ( two - ys ) * ( two - ys )  
             else
                splines_for_z(bz,iz_grid) = -yss * 0.75_double * ( two - ys ) * ( two - ys )
             end if

          ELSE
             splines_for_z(bz,iz_grid) = zero
          END IF
       END DO
    ENDDO

    ncover_yz=BCS_blocks%ncovery*BCS_blocks%ncoverz
    igrid=0


    call start_timer(tmr_std_allocation)
    allocate(send_array(naba_blk_supp%no_naba_blk(iprim)*nsf*n_pts_in_block),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating send_array in do_blip_grad_transform: ",&
         naba_blk_supp%no_naba_blk(iprim)*nsf*n_pts_in_block)
    call reg_alloc_mem(area_basis,naba_blk_supp%no_naba_blk(iprim)*nsf*n_pts_in_block,type_dbl)
    call stop_timer(tmr_std_allocation)
    send_array(:) = 0
    DO naba_blk=1,naba_blk_supp%no_naba_blk(iprim)! naba blks in NOPG order
       ! iprim : primary seq. no. of the atom
       ind_blk=naba_blk_supp%list_naba_blk(naba_blk,iprim) !CC in BCS_blocks
       ind_blk=ind_blk-1

       nx_blk=ind_blk/ncover_yz
       ny_blk=(ind_blk-nx_blk*ncover_yz)/BCS_blocks%ncoverz
       nz_blk=ind_blk-nx_blk*ncover_yz-ny_blk*BCS_blocks%ncoverz

       nx_blk=nx_blk-BCS_blocks%nspanlx+BCS_blocks%nx_origin
       ny_blk=ny_blk-BCS_blocks%nspanly+BCS_blocks%ny_origin
       nz_blk=nz_blk-BCS_blocks%nspanlz+BCS_blocks%nz_origin

       do iz=1,nblkz    ! z-direction in a block
          z=nblkz*(nz_blk-1)+iz-1
          iz_grid=z-nzmin_grid+1
          if(iz_grid > 2*blip_info(spec)%Extent+1) then
             call cq_abort(' ERROR in do_blip_grad_transform_new iz_grid = ', iz_grid,2*blip_info(spec)%Extent+1)
          endif

          imin=imin_for_z(iz_grid)
          imax=imax_for_z(iz_grid)
          do iy=1,nblky    ! y-direction in a block
             y=nblky*(ny_blk-1)+iy-1
             iy_grid=y-nymin_grid+1
             if(iy_grid > 2*blip_info(spec)%Extent+1) then
                call cq_abort(' ERROR in do_blip_grad_transform_new iy_grid = ', iy_grid,2*blip_info(spec)%Extent+1)
             endif

             do ix=1,nblkx    ! x-direction in a block
                x=nblkx*(nx_blk-1)+ix-1
                ix_grid=x-nxmin_grid+1
                if(ix_grid > 2*blip_info(spec)%Extent+1) then
                   call cq_abort(' ERROR in do_blip_grad_transform_new ix_grid = ', ix_grid,2*blip_info(spec)%Extent+1)
                endif
                igrid=igrid+1  ! seq. no. of integration grids
                dsum = zero
                DO bz = imin,imax
                   ys = splines_for_z(bz,iz_grid)
                   do nsf1 = 1,NSF
                      dsum(nsf1)=dsum(nsf1)+ys*inter_2(nsf1,bz,ix_grid,iy_grid)
                   enddo
                END DO
                do nsf1 = 1,NSF
                   send_array(NSF*igrid-NSF+nsf1) = dsum(nsf1)
                enddo
             enddo    ! x-direction in a block
          enddo    ! y-direction in a block
       enddo    ! z-direction in a block
    END DO    ! Loop over naba blocks for this primary atom
    call start_timer(tmr_std_allocation)
    deallocate(inter_1,inter_2,STAT=stat)
    if(stat/=0) call cq_abort('ERROR deallocating memory of inter_1&2 !',stat)
    call reg_dealloc_mem(area_basis,nsf*(2*blip_info(spec)%BlipArraySize+1)*(2*blip_info(spec)%BlipArraySize+1)*&
         (2*blip_info(spec)%Extent+1),type_dbl)
    call reg_dealloc_mem(area_basis,nsf*(2*blip_info(spec)%BlipArraySize+1)*&
         (2*blip_info(spec)%Extent+1)*(2*blip_info(spec)%Extent+1),type_dbl)
    deallocate(splines,splines_for_z,STAT=stat)
    call reg_dealloc_mem(area_basis,(2*blip_info(spec)%BlipArraySize+1)*(2*blip_info(spec)%Extent+2),type_dbl)
    if(stat/=0) call cq_abort('Error deallocating memory to splines !',stat)
    deallocate(imin_for_z,imax_for_z,STAT=stat)
    call reg_dealloc_mem(area_basis,2*(2*blip_info(spec)%Extent+1),type_dbl)
    call stop_timer(tmr_std_allocation)
    return
  end subroutine do_blip_grad_transform_new
!!***

!!****f* blip_grad_transform_module/blip_to_gradgrad_new *
!!
!!  NAME 
!!   blip_to_gradgrad_new
!!  USAGE
!! 
!!  PURPOSE
!!   Does the second derivative of blip coeffients onto the grid
!!  INPUTS
!! 
!! 
!!  USES
!!   do_blip_gradgrad_transform_new
!!   distribute_result
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!!   16/01/2001
!!  MODIFICATION HISTORY
!!   17/05/2001 dave
!!    Added ROBODoc header, indented and stripped subroutine call
!!   11/06/2001 dave
!!    Added GenComms
!!   2006/03/04 06:55 dave
!!    Changed to use new version of distribute result
!!   2008/06/10 ast
!!    Added timers
!!   2011/11/15 17:00 dave
!!    Removed blip module and variables passed to transform
!!  TODO
!!   I think this subroutine can include the work by 
!!   blip_to_grad_transform TM
!!  SOURCE
!!
  subroutine blip_to_gradgrad_new(myid, direction1, direction2, gg_support)

    use datatypes
    use numbers
    use blip, ONLY: blip_info
    use primary_module, ONLY: bundle
    use GenComms, ONLY: my_barrier
    use functions_on_grid, ONLY: gridfunctions
    use species_module, ONLY: nsf_species
    use support_spec_format, ONLY: supports_on_atom

    implicit none

    ! Passed variables
    integer,intent(in) :: myid, direction1,direction2, gg_support

    ! Local variables
    integer :: n_d(3)
    ! n_d(x,y,z)  rank of derivatives with respect to 
    !            x, y and z directions.
    !             At present, sum of n_d should be 2.

    integer :: iprim, nsf_send, spec
    integer :: ii

    call start_timer(tmr_std_basis)
    ! Sets n_d
    n_d=0
    n_d(direction1) = 1
    n_d(direction2) = n_d(direction2) + 1

    gridfunctions(gg_support)%griddata = zero

    do iprim = 1,bundle%mx_iprim
       call my_barrier()
       if( iprim <= bundle%n_prim ) then
          spec = bundle%species(iprim)
          nsf_send = nsf_species(spec)
          call do_blip_gradgrad_transform_new( n_d, iprim, supports_on_atom(iprim), nsf_send)
       else
          nsf_send = 0
       endif
       call my_barrier()
       call distribute_result(myid,iprim,nsf_send,gg_support)
    end do
    call stop_timer(tmr_std_basis)
    return
  end subroutine blip_to_gradgrad_new
!!***

!!****f* blip_grad_transform_module/do_blip_gradgrad_transform_new *
!!
!!  NAME 
!!   do_blip_gradgrad_transform_new
!!  USAGE
!! 
!!  PURPOSE
!!   Does the actual second derivative transform onto a local grid
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!!   16/01/2001
!!  MODIFICATION HISTORY
!!   17/05/2001 dave
!!    Removed NSF=4 dependencies
!!   2008/06/10 ast
!!    Added timers
!!  SOURCE
!!
  subroutine do_blip_gradgrad_transform_new(n_d, iprim, data_blip, nsf)

    use datatypes
    use numbers
    use global_module,   ONLY: rcellx,rcelly,rcellz, area_basis
    use group_module,    ONLY: blocks
    use primary_module,  ONLY: bundle
    use cover_module,    ONLY: BCS_blocks
    use set_blipgrid_module, ONLY: naba_blk_supp
    use comm_array_module,   ONLY: send_array
    use block_module,    ONLY: nx_in_block,ny_in_block,nz_in_block, n_pts_in_block
    use GenComms,        ONLY: cq_abort
    use support_spec_format, ONLY: support_function
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use blip, ONLY: blip_info

    implicit none

    integer,intent(in)     :: n_d(3)
    integer,intent(in)     :: iprim
    integer,intent(in)     :: NSF
    type(support_function),intent(in):: data_blip

    ! Local variables
    real(double),allocatable:: inter_1(:,:,:,:)
    real(double),allocatable:: inter_2(:,:,:,:)
    real(double):: rec_sgs, dx, dy, dz, dxxx, dyyy, dzzz
    integer     :: x, y, z, bx, by, bz 
    integer     :: imin,imax, nsf1, spec

    real(double):: ys, yss, dsum(NSF)

    !NEW VARIABLES introduced by TM
    real(double), allocatable, dimension(:) :: splines
    real(double), allocatable, dimension(:,:) :: splines_for_z
    integer, allocatable, dimension(:) :: imin_for_z, imax_for_z
    integer     :: ierr
    integer     :: ix_grid,iy_grid,iz_grid,ix,iy,iz
    integer     :: nblkx,nblky,nblkz  ! nos. of integ. pts in a block 
    ! along x,y and z direction
    integer     :: naba_blk,nx_blk,ny_blk,nz_blk,igrid
    real(double) :: dcellx_block,dcelly_block,dcellz_block
    real(double) :: dcellx_grid, dcelly_grid, dcellz_grid

    integer     :: ind_blk,ind_blip
    integer     :: nxmin_grid,nxmax_grid,nymin_grid,nymax_grid,&
         nzmin_grid,nzmax_grid
    integer     :: ncover_yz
    integer     :: stat=0,ii

    !real(double) :: sum_send_array1,sum_send_array2,&
         !                sum_send_array3,sum_send_array4
    ! inter_1,inter_2,splines,splines_for_z,imin_for_z,imax_for_z
    ! are dynamically allocated variables.

    spec = bundle%species(iprim)
    rec_sgs = one / blip_info(spec)%SupportGridSpacing
    nblkx=nx_in_block
    nblky=ny_in_block
    nblkz=nz_in_block

    dcellx_block=rcellx/blocks%ngcellx; dcellx_grid=dcellx_block/nblkx
    dcelly_block=rcelly/blocks%ngcelly; dcelly_grid=dcelly_block/nblky
    dcellz_block=rcellz/blocks%ngcellz; dcellz_grid=dcellz_block/nblkz

    nxmin_grid=(naba_blk_supp%nxmin(iprim)-1)*nblkx  
    nxmax_grid= naba_blk_supp%nxmax(iprim)*nblkx-1
    nymin_grid=(naba_blk_supp%nymin(iprim)-1)*nblky  
    nymax_grid= naba_blk_supp%nymax(iprim)*nblky-1
    nzmin_grid=(naba_blk_supp%nzmin(iprim)-1)*nblkz  
    nzmax_grid= naba_blk_supp%nzmax(iprim)*nblkz-1
    !check extent
    ierr=0
    if(2*blip_info(spec)%Extent < nxmax_grid-nxmin_grid) then
       write(io_lun,*) ' ERROR in do_blip_transform : extent for x direction ::', &
            ' extent nxmin_grid nxmax_grid = ',blip_info(spec)%Extent,nxmin_grid,nxmax_grid
       ierr=1
    endif
    if(2*blip_info(spec)%Extent < nymax_grid-nymin_grid) then
       write(io_lun,*) ' ERROR in do_blip_transform : extent for y direction ::', &
            ' extent nymin_grid nymax_grid = ',blip_info(spec)%Extent,nymin_grid,nymax_grid
       ierr=1
    endif
    if(2*blip_info(spec)%Extent < nzmax_grid-nzmin_grid) then
       write(io_lun,*) ' ERROR in do_blip_transform : extent for z direction ::', &
            ' extent nymin_grid nymax_grid = ',blip_info(spec)%Extent,nzmin_grid,nzmax_grid
       ierr=1
    endif
    if(ierr /= 0) call cq_abort('do_blip_gradgrad_transform ',ierr)
    call start_timer(tmr_std_allocation)
    allocate(inter_1(NSF,-blip_info(spec)%BlipArraySize:blip_info(spec)%BlipArraySize, &
         -blip_info(spec)%BlipArraySize:blip_info(spec)%BlipArraySize,2*blip_info(spec)%Extent+1),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to inter_1 !',stat)
    call reg_alloc_mem(area_basis,nsf*(2*blip_info(spec)%BlipArraySize+1)*(2*blip_info(spec)%BlipArraySize+1)*&
         (2*blip_info(spec)%Extent+1),type_dbl)
    call stop_timer(tmr_std_allocation)
    inter_1=zero
    call start_timer(tmr_std_allocation)
    allocate(inter_2(NSF,-blip_info(spec)%BlipArraySize:blip_info(spec)%BlipArraySize,&
         2*blip_info(spec)%Extent+1,2*blip_info(spec)%Extent+1),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to inter_2 !',stat)
    call reg_alloc_mem(area_basis,nsf*(2*blip_info(spec)%BlipArraySize+1)*&
         (2*blip_info(spec)%Extent+1)*(2*blip_info(spec)%Extent+1),type_dbl)
    allocate(splines( -blip_info(spec)%BlipArraySize:blip_info(spec)%BlipArraySize ), &
         splines_for_z(-blip_info(spec)%BlipArraySize:blip_info(spec)%BlipArraySize,1:2*blip_info(spec)%Extent+1))
    call reg_alloc_mem(area_basis,(2*blip_info(spec)%BlipArraySize+1)*(2*blip_info(spec)%Extent+2),type_dbl)
    allocate(imin_for_z(2*blip_info(spec)%Extent+1),imax_for_z(2*blip_info(spec)%Extent+1),STAT=stat)
    call reg_alloc_mem(area_basis,2*(2*blip_info(spec)%Extent+1),type_dbl)
    call stop_timer(tmr_std_allocation)
    inter_2=zero
    imin_for_z = 0
    imax_for_z = 0
    splines = zero
    splines_for_z = zero

    !x,y,z are used for integration grids
    !bx,by,bz are used for support grids
    DO x = nxmin_grid, nxmax_grid
       ix_grid=x-nxmin_grid+1
       dxxx=x*dcellx_grid-bundle%xprim(iprim) ! = dxx+dx in old version
       imin = MAX(-blip_info(spec)%BlipArraySize,int((dxxx-blip_info(spec)%BlipWidth*half)*rec_sgs))
       imax = MIN( blip_info(spec)%BlipArraySize,int((dxxx+blip_info(spec)%BlipWidth*half)*rec_sgs))

       if(imin > blip_info(spec)%BlipArraySize) imin= blip_info(spec)%BlipArraySize
       if(imax <-blip_info(spec)%BlipArraySize) imax=-blip_info(spec)%BlipArraySize

       DO bx = imin, imax
          yss = dxxx - bx * blip_info(spec)%SupportGridSpacing 
          ys = ABS( yss )
          ys = ys * blip_info(spec)%FourOnBlipWidth

          if (yss.ge.zero) then
             yss = blip_info(spec)%FourOnBlipWidth
          else
             yss = -blip_info(spec)%FourOnBlipWidth
          end if
          IF( ys .LE. one ) THEN
             if (n_d(1) == 0) then
                splines( bx ) = one + &
                     ys * ys * ( 0.75_double * ys - 1.5_double )
             elseif(n_d(1) == 1) then
                splines( bx ) = yss * ys * (2.25_double * ys - 3.0_double )
             else
                splines( bx ) = yss * yss * ( 4.5_double * ys - 3.0_double )
             end if
          ELSE IF( ys .LE. two ) THEN
             if (n_d(1) == 0) then
                splines( bx ) = 0.25_double * ( two - ys ) *  &
                     ( two - ys ) * ( two - ys ) 
             elseif(n_d(1) == 1) then
                splines( bx ) = -yss * 0.75_double * ( two - ys ) * ( two - ys )
             else
                splines( bx ) = yss * yss * 1.50_double * ( two - ys )
             end if
          ELSE
             splines( bx ) = zero
          END IF
       END DO
       DO bz = -blip_info(spec)%BlipArraySize, blip_info(spec)%BlipArraySize
          DO by = -blip_info(spec)%BlipArraySize, blip_info(spec)%BlipArraySize
             dsum = zero
             DO bx = imin,imax
                ind_blip = blip_info(spec)%blip_number(bx,by,bz)
                IF( ind_blip .NE. 0 ) THEN
                   ys = splines( bx )
                   do nsf1 = 1,NSF
                      dsum(nsf1) = dsum(nsf1) + ys * data_blip%supp_func(nsf1)%coefficients(ind_blip)
                   enddo
                END IF
             END DO
             do nsf1 = 1,NSF
                inter_1( nsf1,by,bz, ix_grid ) = dsum(nsf1)
             enddo
          END DO
       END DO
    END DO

    DO y = nymin_grid, nymax_grid
       iy_grid=y-nymin_grid+1
       dyyy=y*dcelly_grid-bundle%yprim(iprim)
       imin = MAX(-blip_info(spec)%BlipArraySize,int((dyyy-blip_info(spec)%BlipWidth*half)*rec_sgs))
       imax = MIN( blip_info(spec)%BlipArraySize,int((dyyy+blip_info(spec)%BlipWidth*half)*rec_sgs))

       if(imin > blip_info(spec)%BlipArraySize) imin= blip_info(spec)%BlipArraySize
       if(imax <-blip_info(spec)%BlipArraySize) imax=-blip_info(spec)%BlipArraySize

       DO by = imin,imax
          yss = dyyy - by * blip_info(spec)%SupportGridSpacing 
          ys = ABS( yss )
          ys = ys * blip_info(spec)%FourOnBlipWidth

          if (yss.ge.zero) then
             yss = blip_info(spec)%FourOnBlipWidth
          else
             yss = -blip_info(spec)%FourOnBlipWidth
          end if

          IF( ys .LE. one ) THEN
             if (n_d(2) == 0) then
                splines( by ) = one + ys * ys * ( 0.75_double * ys - 1.5_double )
             else if (n_d(2) == 1) then
                splines( by ) = yss * ys * (2.25_double * ys - 3.0_double )
             else
                splines( by ) = yss * yss * (4.5_double * ys - 3.0_double )
             end if
          ELSE IF( ys .LE. two ) THEN
             if (n_d(2) == 0) then
                splines( by ) = 0.25_double * ( two - ys ) * &
                     ( two - ys ) * ( two - ys )
             else if (n_d(2) == 1) then
                splines( by ) = -yss * 0.75_double * ( two - ys ) * ( two - ys )
             else
                splines( by ) = yss * yss * 1.5_double * ( two - ys )
             end if
          ELSE
             splines( by ) = zero
          END IF
       END DO
       DO x = nxmin_grid, nxmax_grid
          ix_grid=x-nxmin_grid+1
          DO bz = -blip_info(spec)%BlipArraySize, blip_info(spec)%BlipArraySize
             dsum = zero
             DO by = imin,imax
                ys = splines( by )
                do nsf1 = 1,NSF
                   dsum(nsf1) = dsum(nsf1) + ys*inter_1( nsf1,by,bz,ix_grid)
                enddo
             END DO
             do nsf1 = 1,NSF
                inter_2( nsf1, bz,ix_grid,iy_grid) = dsum(nsf1)
             enddo
          END DO
       END DO
    END DO

    DO z = nzmin_grid, nzmax_grid
       iz_grid=z-nzmin_grid+1
       dzzz=z*dcellz_grid-bundle%zprim(iprim)
       imin_for_z(iz_grid) = MAX(-blip_info(spec)%BlipArraySize,&
            int((dzzz-blip_info(spec)%BlipWidth*half)*rec_sgs))
       imax_for_z(iz_grid) = MIN( blip_info(spec)%BlipArraySize,&
            int((dzzz+blip_info(spec)%BlipWidth*half)*rec_sgs))

       if(imin_for_z(iz_grid) > blip_info(spec)%BlipArraySize) &
            imin_for_z(iz_grid) = blip_info(spec)%BlipArraySize
       if(imax_for_z(iz_grid) <-blip_info(spec)%BlipArraySize) &
            imax_for_z(iz_grid) =-blip_info(spec)%BlipArraySize

       DO bz = imin_for_z(iz_grid),imax_for_z(iz_grid)
          yss = dzzz - bz * blip_info(spec)%SupportGridSpacing 
          ys = ABS( yss )
          ys = ys * blip_info(spec)%FourOnBlipWidth
          if (yss.ge.zero) then
             yss = blip_info(spec)%FourOnBlipWidth
          else
             yss =-blip_info(spec)%FourOnBlipWidth
          end if
          IF( ys .LE. one ) THEN
             if (n_d(3) == 0) then
                splines_for_z(bz,iz_grid) = one +                      &
                     ys * ys * ( 0.75_double * ys - 1.5_double )
             elseif (n_d(3) == 1) then
                splines_for_z(bz,iz_grid) = yss * ys * (2.25_double * ys - 3.0_double )
             else
                splines_for_z(bz,iz_grid) = yss * yss * (4.5_double * ys - 3.0_double )
             end if
          ELSE IF( ys .LE. two ) THEN
             if (n_d(3) == 0) then
                splines_for_z(bz,iz_grid) = 0.25_double * ( two - ys ) *    &
                     ( two - ys ) * ( two - ys )  
             elseif (n_d(3) == 1) then
                splines_for_z(bz,iz_grid) = -yss * 0.75_double * &
                     ( two - ys ) * ( two - ys )
             else
                splines_for_z(bz,iz_grid) =  yss * yss * 1.5_double * ( two - ys )
             end if
          ELSE
             splines_for_z(bz,iz_grid) = zero
          END IF
       END DO
    ENDDO

    ncover_yz=BCS_blocks%ncovery*BCS_blocks%ncoverz
    igrid=0

 
    call start_timer(tmr_std_allocation)
    allocate(send_array(naba_blk_supp%no_naba_blk(iprim)*nsf*n_pts_in_block),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating send_array in do_blip_gradgrad_transform: ",&
         naba_blk_supp%no_naba_blk(iprim)*nsf*n_pts_in_block)
    call reg_alloc_mem(area_basis,naba_blk_supp%no_naba_blk(iprim)*nsf*n_pts_in_block,type_dbl)
    call stop_timer(tmr_std_allocation)
    send_array(:) = zero
    DO naba_blk=1,naba_blk_supp%no_naba_blk(iprim)! naba blks in NOPG order
       ! iprim : primary seq. no. of the atom
       ind_blk=naba_blk_supp%list_naba_blk(naba_blk,iprim) !CC in BCS_blocks
       ind_blk=ind_blk-1

       nx_blk=ind_blk/ncover_yz
       ny_blk=(ind_blk-nx_blk*ncover_yz)/BCS_blocks%ncoverz
       nz_blk=ind_blk-nx_blk*ncover_yz-ny_blk*BCS_blocks%ncoverz

       nx_blk=nx_blk-BCS_blocks%nspanlx+BCS_blocks%nx_origin
       ny_blk=ny_blk-BCS_blocks%nspanly+BCS_blocks%ny_origin
       nz_blk=nz_blk-BCS_blocks%nspanlz+BCS_blocks%nz_origin

       do iz=1,nblkz    ! z-direction in a block
          z=nblkz*(nz_blk-1)+iz-1
          iz_grid=z-nzmin_grid+1
          if(iz_grid > 2*blip_info(spec)%Extent+1) then
             call cq_abort(' ERROR in do_blip_gradgrad_transform_new iz_grid = ', iz_grid,2*blip_info(spec)%Extent+1)
          endif

          imin=imin_for_z(iz_grid)
          imax=imax_for_z(iz_grid)
          do iy=1,nblky    ! y-direction in a block
             y=nblky*(ny_blk-1)+iy-1
             iy_grid=y-nymin_grid+1
             if(iy_grid > 2*blip_info(spec)%Extent+1) then
                call cq_abort(' ERROR in do_blip_gradgrad_transform_new iy_grid = ', iy_grid,2*blip_info(spec)%Extent+1)
             endif

             do ix=1,nblkx    ! x-direction in a block
                x=nblkx*(nx_blk-1)+ix-1
                ix_grid=x-nxmin_grid+1
                if(ix_grid > 2*blip_info(spec)%Extent+1) then
                   call cq_abort(' ERROR in do_blip_gradgrad_transform_new ix_grid = ', ix_grid,2*blip_info(spec)%Extent+1)
                endif

                igrid=igrid+1  ! seq. no. of integration grids
                dsum = zero
                DO bz = imin,imax
                   ys = splines_for_z(bz,iz_grid)
                   do nsf1 = 1,NSF
                      dsum(nsf1)=dsum(nsf1)+ys*inter_2(nsf1,bz,ix_grid,iy_grid)
                   enddo
                END DO
                do nsf1 = 1,NSF
                   send_array(NSF*igrid-NSF+nsf1) = dsum(nsf1)
                enddo
             enddo    ! x-direction in a block
          enddo    ! y-direction in a block
       enddo    ! z-direction in a block
    END DO    ! Loop over naba blocks for this primary atom
    call start_timer(tmr_std_allocation)
    deallocate(inter_1,inter_2,STAT=stat)
    if(stat/=0) call cq_abort(' ERROR in do_blip_gradgrad_transform_new: deallocation')
    call reg_dealloc_mem(area_basis,nsf*(2*blip_info(spec)%BlipArraySize+1)*(2*blip_info(spec)%BlipArraySize+1)*&
         (2*blip_info(spec)%Extent+1),type_dbl)
    call reg_dealloc_mem(area_basis,nsf*(2*blip_info(spec)%BlipArraySize+1)*&
         (2*blip_info(spec)%Extent+1)*(2*blip_info(spec)%Extent+1),type_dbl)
    deallocate(splines,splines_for_z,STAT=stat)
    call reg_dealloc_mem(area_basis,(2*blip_info(spec)%BlipArraySize+1)*(2*blip_info(spec)%Extent+2),type_dbl)
    if(stat/=0) call cq_abort('Error deallocating memory to splines !',stat)
    deallocate(imin_for_z,imax_for_z,STAT=stat)
    call reg_dealloc_mem(area_basis,2*(2*blip_info(spec)%Extent+1),type_dbl)
    call stop_timer(tmr_std_allocation)
    return
  end subroutine do_blip_gradgrad_transform_new
!!***

!!****f* blip_grid_transform_module/inverse_blip_transform_new *
!!
!!  NAME 
!!   inverse_blip_transform_new
!!  USAGE
!! 
!!  PURPOSE
!!   Takes the support functions on the grid and returns the 
!!   corresponding blip coefficients
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!!   27/06/2000
!!  MODIFICATION HISTORY
!!   17/05/2001 dave
!!    Added ROBODoc header, stripped variables list
!!   11/06/2001 dave
!!    Added GenComms
!!   2008/06/10 ast
!!    Added timers
!!   2011/11/15 17:00 dave
!!    Removed blip module and variables passed to transform
!!  SOURCE
!!
  subroutine inverse_blip_transform_new(myid,dsupport, data_dblip, n_prim)

    use datatypes
    use numbers
    use primary_module, ONLY:bundle
    use comm_array_module, ONLY: send_array
    use set_blipgrid_module, ONLY:naba_blk_supp,comBG
    use block_module,        ONLY: n_pts_in_block  ! = blocks%nm_group(:)
    use blip, ONLY: blip_info
    use GenComms, ONLY: my_barrier
    use species_module, ONLY: nsf_species
    use support_spec_format, ONLY: support_function

    implicit none

    integer,intent(in) :: myid, dsupport, n_prim
    type(support_function), dimension(n_prim) :: data_dblip

    ! local variables
    integer :: iprim, nsf_recv, spec

    call start_timer(tmr_std_basis)
    do iprim = 1, bundle%mx_iprim
       !  dsupport shows values on integration grids
       !  Bundle responsible nodes accumulates their 
       ! responsible values on integration grids from
       ! Domain responsible nodes and put them into 
       ! send_array. ( this array is defined in comm_array_module.)

       if(iprim<=bundle%n_prim) then 
          spec = bundle%species(iprim)
          nsf_recv = nsf_species(spec)
       else
          nsf_recv = 0
       end if
       !send_array(:)=zero
       call collect_result (myid,iprim,nsf_recv,dsupport)

       call my_barrier()   ! this is not needed.

       if(iprim <= bundle%n_prim) then
          call do_inverse_blip_new ( myid, iprim, data_dblip(iprim), nsf_recv)
       endif
    end do
    call stop_timer(tmr_std_basis)
    return
  end subroutine inverse_blip_transform_new
!!***

!!****f* blip_grid_transform_module/do_inverse_blip_new *
!!
!!  NAME 
!!   do_inverse_blip_new
!!  USAGE
!! 
!!  PURPOSE
!!   Performs the transform; new version:
!!   1. to be consistent with new set_grid
!!   2. the concept of primary set and covering set is introduced.
!!   3. the present version assumes NSF is 4 for all atoms
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!!   23/06/2000
!!  MODIFICATION HISTORY
!!   17/05/2001 dave
!!    Removed NSF=4 dependencies
!!   2008/06/10 ast
!!    Added timers
!!  SOURCE
!!
  subroutine do_inverse_blip_new( myid, iprim, data_dblip, nsf)

    ! modules
    use datatypes
    use numbers
    use global_module,       ONLY:rcellx,rcelly,rcellz, area_basis
    use group_module,        ONLY:blocks
    use primary_module,      ONLY:bundle
    use cover_module,        ONLY:BCS_blocks
    use set_blipgrid_module, ONLY:comBG,naba_blk_supp
    use comm_array_module,   ONLY:send_array

    use block_module,ONLY: n_pts_in_block,nx_in_block,ny_in_block,nz_in_block
    use GenComms, ONLY: cq_abort
    use support_spec_format, ONLY: support_function
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use blip, ONLY: blip_info

    ! Shared variable
    implicit none

    integer,intent(in):: iprim,myid
    integer,intent(in):: nsf
    type(support_function),intent(out):: data_dblip

    ! Local variables
    real(double),allocatable:: inter_1(:,:,:,:)
    real(double),allocatable:: inter_2(:,:,:,:)
    real(double):: rec_sgs, dxxx, dyyy, dzzz
    integer     :: x, y, z, bx, by, bz, imin, imax, nsf1, spec
    integer     :: nxmin_grid,nxmax_grid
    integer     :: nymin_grid,nymax_grid
    integer     :: nzmin_grid,nzmax_grid
    integer, allocatable, dimension(:) :: imin_for_z, imax_for_z

    real(double), allocatable, dimension(:) :: splines
    real(double), allocatable, dimension(:,:) :: splines_for_z
    real(double):: ys, yss, dsum(NSF)

    integer     :: ierr,stat
    integer     :: ix_grid,iy_grid,iz_grid,ix,iy,iz
    real(double):: dcellx_block,dcelly_block,dcellz_block
    real(double):: dcellx_grid, dcelly_grid, dcellz_grid
    integer     :: ind_blip,igrid,ind_blk,nx_blk,ny_blk,nz_blk
    integer     :: ncover_yz,naba_blk

    !!   nblkx,nblky,nblkz : number of integration grid points in a block
    !!     this variables should be passed from (blocks)
    integer :: nblkx,nblky,nblkz

    ! For checking
    !      integer   :: isf,ii
    !      real(double) :: sum
    !**************************************************************************
    !     Start of subroutine

    spec = bundle%species(iprim)
    call start_timer(tmr_std_allocation)
    allocate(inter_1(NSF,-blip_info(spec)%BlipArraySize:blip_info(spec)%BlipArraySize,&
         -blip_info(spec)%BlipArraySize:blip_info(spec)%BlipArraySize,2*blip_info(spec)%Extent+1),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to inter_1 !',stat)
    call reg_alloc_mem(area_basis,nsf*(2*blip_info(spec)%BlipArraySize+1)*(2*blip_info(spec)%BlipArraySize+1)*&
         (2*blip_info(spec)%Extent+1),type_dbl)
    call stop_timer(tmr_std_allocation)
    inter_1=zero
    call start_timer(tmr_std_allocation)
    allocate(inter_2(NSF,-blip_info(spec)%BlipArraySize:blip_info(spec)%BlipArraySize,&
         2*blip_info(spec)%Extent+1 ,2*blip_info(spec)%Extent+1),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to inter_2 !',stat)
    call reg_alloc_mem(area_basis,nsf*(2*blip_info(spec)%BlipArraySize+1)*&
         (2*blip_info(spec)%Extent+1)*(2*blip_info(spec)%Extent+1),type_dbl)
    allocate(splines( -blip_info(spec)%BlipArraySize:blip_info(spec)%BlipArraySize ), &
         splines_for_z(-blip_info(spec)%BlipArraySize:blip_info(spec)%BlipArraySize,1:2*blip_info(spec)%Extent+1))
    call reg_alloc_mem(area_basis,(2*blip_info(spec)%BlipArraySize+1)*(2*blip_info(spec)%Extent+2),type_dbl)
    allocate(imin_for_z(2*blip_info(spec)%Extent+1),imax_for_z(2*blip_info(spec)%Extent+1),STAT=stat)
    call reg_alloc_mem(area_basis,2*(2*blip_info(spec)%Extent+1),type_dbl)
    call stop_timer(tmr_std_allocation)
    inter_2=zero
    imin_for_z = 0
    imax_for_z = 0
    splines = zero
    splines_for_z = zero

    rec_sgs = one / blip_info(spec)%SupportGridSpacing
    nblkx=nx_in_block
    nblky=ny_in_block
    nblkz=nz_in_block

    dcellx_block=rcellx/blocks%ngcellx; dcellx_grid=dcellx_block/nblkx
    dcelly_block=rcelly/blocks%ngcelly; dcelly_grid=dcelly_block/nblky
    dcellz_block=rcellz/blocks%ngcellz; dcellz_grid=dcellz_block/nblkz

    nxmin_grid=(naba_blk_supp%nxmin(iprim)-1)*nblkx
    nxmax_grid= naba_blk_supp%nxmax(iprim)*nblkx-1
    nymin_grid=(naba_blk_supp%nymin(iprim)-1)*nblky
    nymax_grid= naba_blk_supp%nymax(iprim)*nblky-1
    nzmin_grid=(naba_blk_supp%nzmin(iprim)-1)*nblkz
    nzmax_grid= naba_blk_supp%nzmax(iprim)*nblkz-1

    !check extent
    ierr=0
    if(2*blip_info(spec)%Extent < nxmax_grid-nxmin_grid) then
       write(io_lun,*) myid,' ERROR in do_inverse_blip : extent for x dir ::', &
            ' extent nxmin_grid nxmax_grid = ',blip_info(spec)%Extent,nxmin_grid,nxmax_grid
       ierr=1
    endif
    if(2*blip_info(spec)%Extent < nymax_grid-nymin_grid) then
       write(io_lun,*) myid,' ERROR in do_inverse_blip : extent for y dir ::', &
            ' extent nymin_grid nymax_grid = ',blip_info(spec)%Extent,nymin_grid,nymax_grid
       ierr=1
    endif
    if(2*blip_info(spec)%Extent < nzmax_grid-nzmin_grid) then
       write(io_lun,*) myid,' ERROR in do_inverse_blip : extent for z dir ::', &
            ' extent nymin_grid nymax_grid = ',blip_info(spec)%Extent,nzmin_grid,nzmax_grid
       ierr=1
    endif
    if(ierr /= 0) call cq_abort('do_inverse_blip ',ierr)
    DO z = nzmin_grid , nzmax_grid
       iz_grid=z-nzmin_grid+1
       dzzz=z*dcellz_grid-bundle%zprim(iprim)
       imin_for_z(iz_grid)=MAX(-blip_info(spec)%BlipArraySize, &
            int((dzzz-blip_info(spec)%BlipWidth*half)*rec_sgs))
       imax_for_z(iz_grid)=MIN( blip_info(spec)%BlipArraySize, &
            int((dzzz+blip_info(spec)%BlipWidth*half)*rec_sgs))
       if(imin_for_z(iz_grid) > blip_info(spec)%BlipArraySize) then
          imin_for_z(iz_grid)=blip_info(spec)%BlipArraySize
       endif
       if(imax_for_z(iz_grid) < -blip_info(spec)%BlipArraySize) then
          imax_for_z(iz_grid)=-blip_info(spec)%BlipArraySize
       endif

       DO bz = imin_for_z(iz_grid), imax_for_z(iz_grid)
          yss = dzzz - bz * blip_info(spec)%SupportGridSpacing 
          ys = ABS( yss )
          ys = ys * blip_info(spec)%FourOnBlipWidth
          IF( ys .LE. one ) THEN
             splines_for_z(bz,iz_grid) = one +                     &
                  ys * ys * ( 0.75_double * ys - 1.5_double )
          ELSE IF( ys .LE. two ) THEN
             splines_for_z(bz,iz_grid) = 0.25_double * ( two - ys ) *   &
                  ( two - ys ) * ( two - ys )  
          ELSE
             splines_for_z(bz,iz_grid) = zero
          END IF
       END DO
    ENDDO

    ncover_yz=BCS_blocks%ncovery*BCS_blocks%ncoverz
    igrid=0
    DO naba_blk=1,naba_blk_supp%no_naba_blk(iprim)! naba blks (in NOPG order)
       ! iprim : prim. seq. # of atm
       ind_blk=naba_blk_supp%list_naba_blk(naba_blk,iprim) !CC in BCS_blocks
       ind_blk=ind_blk-1

       nx_blk=ind_blk/ncover_yz
       ny_blk=(ind_blk-nx_blk*ncover_yz)/BCS_blocks%ncoverz
       nz_blk=ind_blk-nx_blk*ncover_yz-ny_blk*BCS_blocks%ncoverz

       nx_blk=nx_blk-BCS_blocks%nspanlx+BCS_blocks%nx_origin
       ny_blk=ny_blk-BCS_blocks%nspanly+BCS_blocks%ny_origin
       nz_blk=nz_blk-BCS_blocks%nspanlz+BCS_blocks%nz_origin

       do iz=1,nblkz    ! z-direction in a block
          z=nblkz*(nz_blk-1)+iz-1
          iz_grid=z-nzmin_grid+1

          imin=imin_for_z(iz_grid)
          imax=imax_for_z(iz_grid)
          do iy=1,nblky    ! y-direction in a block
             y=nblky*(ny_blk-1)+iy-1
             iy_grid=y-nymin_grid+1
             if(iy_grid < 1 .OR. iy_grid > nymax_grid-nymin_grid+1) then
                call cq_abort(' ERROR iy_grid = ',iy_grid,nymax_grid-nymin_grid+1)
             endif
             do ix=1,nblkx    ! x-direction in a block
                x=nblkx*(nx_blk-1)+ix-1
                ix_grid=x-nxmin_grid+1
                if(ix_grid < 1 .OR. ix_grid > nxmax_grid-nxmin_grid+1) then
                   call cq_abort(' ERROR ix_grid = ',ix_grid,nxmax_grid-nxmin_grid+1)
                endif

                igrid=igrid+1  ! seq. no. of integration grids

                if(igrid > naba_blk_supp%no_naba_blk(iprim)*n_pts_in_block .OR. igrid < 1) then
                   call cq_abort('ERROR!!!!!!!   iprim,naba_blk,iz,iy,ix = ',&
                        igrid,naba_blk_supp%no_naba_blk(iprim)*n_pts_in_block)
                endif
                do nsf1 = 1,NSF  ! data_naba_grid(1,igrid)
                   dsum(nsf1) = send_array(NSF*igrid-NSF+nsf1)
                enddo
                DO bz = imin,imax
                   ys = splines_for_z(bz,iz_grid)
                   do nsf1 = 1,NSF
                      inter_2(nsf1,bz,ix_grid,iy_grid)= &
                           inter_2( nsf1,bz,ix_grid,iy_grid)+dsum(nsf1)*ys
                   enddo
                END DO
             enddo    ! x-direction in a block
          enddo    ! y-direction in a block
       enddo    ! z-direction in a block
    END DO    ! Loop over naba blocks for this primary atom
    DO y = nymin_grid, nymax_grid
       iy_grid=y-nymin_grid+1
       dyyy=y*dcelly_grid-bundle%yprim(iprim)
       imin = MAX(-blip_info(spec)%BlipArraySize,int((dyyy-blip_info(spec)%BlipWidth*half)*rec_sgs))
       imax = MIN( blip_info(spec)%BlipArraySize,int((dyyy+blip_info(spec)%BlipWidth*half)*rec_sgs))
       if(imin >  blip_info(spec)%BlipArraySize) imin= blip_info(spec)%BlipArraySize
       if(imax < -blip_info(spec)%BlipArraySize) imax=-blip_info(spec)%BlipArraySize
       DO by = imin, imax
          yss = dyyy - by * blip_info(spec)%SupportGridSpacing 
          ys = ABS( yss )
          ys = ys * blip_info(spec)%FourOnBlipWidth
          IF( ys .LE. one ) THEN
             splines( by ) = one +  &
                  ys * ys * ( 0.75_double * ys - 1.5_double )
          ELSE IF( ys .LE. two ) THEN
             splines( by ) = 0.25_double * ( two - ys ) *  &
                  ( two - ys ) * ( two - ys )  
          ELSE
             splines( by ) = zero
          END IF
       END DO
       DO x = nxmin_grid, nxmax_grid
          ix_grid=x-nxmin_grid+1
          DO bz = -blip_info(spec)%BlipArraySize, blip_info(spec)%BlipArraySize
             do nsf1 = 1,NSF
                dsum(nsf1) = inter_2(nsf1,bz,ix_grid,iy_grid )
             enddo
             DO by = imin, imax
                ys = splines( by )
                do nsf1 = 1,NSF
                   inter_1(nsf1,by,bz,ix_grid)=dsum(nsf1)*ys + &
                        inter_1(nsf1,by,bz,ix_grid)
                enddo
             END DO
          END DO
       END DO
    END DO
    DO x = nxmin_grid, nxmax_grid
       ix_grid=x-nxmin_grid+1
       dxxx=x*dcellx_grid-bundle%xprim(iprim)
       imin = MAX(-blip_info(spec)%BlipArraySize,int((dxxx-blip_info(spec)%BlipWidth*half)*rec_sgs))
       imax = MIN( blip_info(spec)%BlipArraySize,int((dxxx+blip_info(spec)%BlipWidth*half)*rec_sgs))
       if(imin >  blip_info(spec)%BlipArraySize) imin= blip_info(spec)%BlipArraySize
       if(imax < -blip_info(spec)%BlipArraySize) imax=-blip_info(spec)%BlipArraySize
       DO bx =imin, imax
          yss = dxxx - bx * blip_info(spec)%SupportGridSpacing
          ys = ABS( yss )
          ys = ys * blip_info(spec)%FourOnBlipWidth
          IF( ys .LE. one ) THEN
             splines( bx ) = one +  &
                  ys * ys * ( 0.75_double * ys - 1.5_double )
          ELSE IF( ys .LE. two ) THEN
             splines( bx ) = 0.25_double * ( two - ys ) *  &
                  ( two - ys ) * ( two - ys ) 
          ELSE
             splines( bx ) = zero
          END IF
       END DO
       DO bz = -blip_info(spec)%BlipArraySize, blip_info(spec)%BlipArraySize
          DO by = -blip_info(spec)%BlipArraySize, blip_info(spec)%BlipArraySize
             do nsf1 = 1,NSF
                dsum(nsf1) = inter_1( nsf1,by,bz,ix_grid)
             enddo
             DO bx = imin, imax
                ind_blip = blip_info(spec)%blip_number(bx,by,bz)
                IF(ind_blip /= 0 ) THEN
                   ys = splines( bx )
                   do nsf1 = 1,NSF
                      data_dblip%supp_func(nsf1)%coefficients(ind_blip)=dsum(nsf1)*ys+ &
                           data_dblip%supp_func(nsf1)%coefficients(ind_blip)
                   enddo
                END IF
             END DO
          END DO
       END DO
    END DO
    call start_timer(tmr_std_allocation)
    deallocate(inter_1,inter_2,STAT=stat)
    if(stat/=0) call cq_abort('Error deallocating memory to inter_1 !',stat)
    call reg_dealloc_mem(area_basis,nsf*(2*blip_info(spec)%BlipArraySize+1)*(2*blip_info(spec)%BlipArraySize+1)*&
         (2*blip_info(spec)%Extent+1),type_dbl)
    call reg_dealloc_mem(area_basis,nsf*(2*blip_info(spec)%BlipArraySize+1)*&
         (2*blip_info(spec)%Extent+1)*(2*blip_info(spec)%Extent+1),type_dbl)
    call reg_dealloc_mem(area_basis,size(send_array),type_dbl)
    deallocate(splines,splines_for_z,STAT=stat)
    call reg_dealloc_mem(area_basis,(2*blip_info(spec)%BlipArraySize+1)*(2*blip_info(spec)%Extent+2),type_dbl)
    if(stat/=0) call cq_abort('Error deallocating memory to splines !',stat)
    deallocate(send_array,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating send_array in do_inverse_blip: ",stat)
    deallocate(imin_for_z,imax_for_z,STAT=stat)
    call reg_dealloc_mem(area_basis,2*(2*blip_info(spec)%Extent+1),type_dbl)
    call stop_timer(tmr_std_allocation)
    return
  end subroutine do_inverse_blip_new
!!***

  subroutine collect_result(myid,iprim,nsf_recv,support)
    !  made from distribute_result 
    !     27/6/2000 Tsuyoshi Miyazaki
    !  collects the values on integration grids to do inverse
    !  blip-grid transforms.
    !
    !  Maybe...
    !  we can use this subroutine even when we treat various nsf.
    !   26/6/2000 Now we assume that this subroutine is called once to 
    !             send blip-grid transformed support functions of one
    !             primary atom. If the number of support func. is
    !             different for different atoms, nunit and nsize varies
    !             at each call.
    !
    !    5/7/2000  I have changed the whole structure of this subroutine.
    !             This subroutine will do the inverse of distribute_result.
    !             Thus, I made this subroutine from distribute_result 
    !             by replacing (MPI_issend + MPI_recv) with 
    !             (MPI_irecv+MPI_ssend+MPI_wait).
    !             I think this change is simple and the present code is 
    !             easy to understand if you can understand distribute_result,
    !             although the name of (send_array) and (recv_array) might
    !             look curious.
    !              As for the memory requirement, the present version is better.
    !             send_array, which is used as the array to store the sent data,
    !             and recv_array, which is used as the one to be sent, have the 
    !             same size as those in blip-grid transforms. Thus, we can reuse
    !             the arrays used in blip-grid transforms.
    !  2006/06/14
    !  My previous expectation in the above is too optimistic....
    !  As we now treat various NSF, nsf for sending and nsf for receiving 
    !  are generally different.

    use datatypes
    use numbers
    use mpi
    use set_blipgrid_module, ONLY: comBG,naba_atm
    use comm_array_module,   ONLY: send_array,recv_array
    use block_module,        ONLY: n_pts_in_block  != blocks%nm_group(:)
    use GenComms, ONLY: my_barrier, cq_abort
    use functions_on_grid, ONLY: gridfunctions, fn_on_grid
    use species_module, ONLY: nsf_species
    use global_module, ONLY: species_glob, sf, area_basis
    use atoms, ONLY: atoms_on_node
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl, type_int

    implicit none
    integer,intent(in) :: myid,iprim,nsf_recv
    integer :: support, msize
    integer :: mynode,nunit_recv, nsf_send, nunit_send,ind_alp_i_blk
    integer :: inode,jnode,nnd_rem,ibegin,nsize,my_ibegin
    integer :: ipair,igrid,isf,prim_blk,naba_atm_tmp,loc1,loc2,ind1,ind2
    real(double),pointer :: recv_ptr(:)

    integer :: tag,ierr,irc,stat
    integer, allocatable, dimension(:) :: nrecv_req
    integer :: nwait_stat(MPI_STATUS_SIZE)
    integer :: isend,off


    mynode=myid+1
    tag=1
    ! 
    !Prepares receiving the values on integration grids
    !(those values will be put into send_array)
    !
    !  if nnd_rem is my node, the following just keeps the first address 
    !  in send_array. This will be pointed to recv_ptr.
    !  

    nsize = 0
    if(comBG%no_recv_node(iprim) > 0) then
       nunit_recv=nsf_recv*n_pts_in_block   
       do inode=1,comBG%no_recv_node(iprim)
          nsize = nsize+comBG%no_naba_blk(inode,iprim)
       enddo
       !write(io_lun,*) 'Allocating send_array: ',nsize
       call start_timer(tmr_std_allocation)
       allocate(send_array(nsize*nunit_recv),STAT=stat)
       if(stat/=0) call cq_abort("Error allocating send_array in collect_result: ", nsize,nunit_recv)
       call reg_alloc_mem(area_basis,nsize*nunit_recv,type_dbl)
       call stop_timer(tmr_std_allocation)
       send_array(:)=zero
       call start_timer(tmr_std_allocation)
       allocate(nrecv_req(comBG%no_recv_node(iprim)), STAT=ierr)
       if(ierr/=0) call cq_abort("Error allocating nrecv_req in collect_result: ",comBG%no_recv_node(iprim),ierr)
       call reg_alloc_mem(area_basis,comBG%no_recv_node(iprim),type_int)
       call stop_timer(tmr_std_allocation)
       do inode=1,comBG%no_recv_node(iprim)
          nnd_rem=comBG%list_recv_node(inode,iprim)
          if(inode == 1) then
             ibegin=1
          else
             ibegin=ibegin+comBG%no_naba_blk(inode-1,iprim)*nunit_recv
          endif

          nsize=comBG%no_naba_blk(inode,iprim)*nunit_recv
          if(nnd_rem == mynode) then
             my_ibegin=ibegin
          else
             call MPI_irecv(send_array(ibegin),nsize,MPI_DOUBLE_PRECISION,&
                  nnd_rem-1,tag,MPI_COMM_WORLD,nrecv_req(inode),ierr)
          endif
       enddo
    endif ! if there are receiving nodes

    call my_barrier() ! this is not needed.

    !Makes Arrays and Sends
    if(comBG%no_send_node(iprim) > 0) then
       nsize=0
       msize=0
       do jnode=1,comBG%no_send_node(iprim)
          nnd_rem=comBG%list_send_node(jnode,iprim)
          nsf_send = nsf_species(species_glob(atoms_on_node(iprim,nnd_rem)))
          nunit_send = nsf_send*n_pts_in_block
          msize=max(msize,comBG%no_sent_pairs(jnode,iprim)*nunit_send)
       end do
       !write(io_lun,*) 'Allocating recv_array: ',nsize
       call start_timer(tmr_std_allocation)
       allocate(recv_array(msize),STAT=stat)
       if(stat/=0) call cq_abort("Error allocating recv_array in collect_result: ",nsize)
       call reg_alloc_mem(area_basis,msize,type_dbl)
       call stop_timer(tmr_std_allocation)
       isend=comBG%ibeg_recv_call(iprim)-1
       off = 0
       do jnode=1,comBG%no_send_node(iprim)
          isend=isend+1
          nnd_rem=comBG%list_send_node(jnode,iprim)
          nsf_send = nsf_species(species_glob(atoms_on_node(iprim,nnd_rem)))
          nunit_send = nsf_send*n_pts_in_block
          nsize=comBG%no_sent_pairs(jnode,iprim)*nunit_send
          if(nnd_rem == mynode) then
             recv_ptr => send_array(my_ibegin:my_ibegin+nsize-1)
          else
             recv_ptr => recv_array
          endif

          do ipair=1,comBG%no_sent_pairs(jnode,iprim) !naba blks of the iprim-th atm
             prim_blk=comBG%table_blk(ipair,isend)     ! primary block in my domain
             naba_atm_tmp=comBG%table_pair(ipair,isend)! naba atm of the primary blk
             loc1= nunit_send*(ipair-1)
             ind_alp_i_blk = naba_atm(sf)%ibegin_blk_orb(prim_blk)+ &
                  naba_atm(sf)%ibeg_orb_atom(naba_atm_tmp, prim_blk) -1 
             ! naba_atm(sf)%ibegin_blk_orb(iprim_blk): initial position of support for the present primary block.
             ! naba_atm(sf)%ibeg_orb_atom(naba, iprim_blk)
             ! : initial position of the present naba atom (naba) in the neighbour-atom orbitals of the 
             ! present primary block.
             !loc2 = n_pts_in_block * (ind_alp_i_blk-1) + 1
             loc2 = n_pts_in_block * (ind_alp_i_blk-1)

             if(loc1+nunit_send>nsize .or. loc1<0) then
                call cq_abort(' ERROR loc1 in collect_result ',nunit_send+loc1, nsize)
             endif
             if(loc2+nunit_send>gridfunctions(support)%size .or. loc2<0) then
                write(io_lun,*) ' ERROR loc2 in collect_result for mynode= ',&
                     mynode,' loc2 = ',loc2+nunit_send,gridfunctions(support)%size
                write(io_lun,*) '  ERROR prim_blk, naba_atm_tmp = ',prim_blk,naba_atm_tmp,&
                     ' ibegin_blk = ',naba_atm(sf)%ibegin_blk(prim_blk),&
                     ' naba_atom_of_blk = ',naba_atm(sf)%no_of_atom(prim_blk)
                call cq_abort("Stopping in inverse_blip_transform")
             endif

             !!------ TRANSPOSE
             do isf=1,nsf_send
                ind2=(isf-1)*n_pts_in_block+loc2
                do igrid=1,n_pts_in_block
                   ind1=(igrid-1)*nsf_send+isf
                   ind2=ind2+1
                   recv_ptr(ind1+loc1)=gridfunctions(support)%griddata(ind2)
                enddo
             enddo
             !  ind1=loc1
             ! do igrid=1,n_pts_in_block
             !  do isf=1,NSF
             !    ind1=ind1+1  !  <=>  ind1=(igrid-1)*NSF+isf + loc1
             !    ind2=(isf-1)*n_pts_in_block+igrid
             !    recv_ptr(ind1)=support(loc2+ind2)
             !  enddo
             ! enddo
             !!------ TRANSPOSE
          enddo ! loop over pairs of block and atoms

          if(nnd_rem /= mynode) then
             call MPI_ssend(recv_array,nsize,MPI_DOUBLE_PRECISION,nnd_rem-1, &
                  tag,MPI_COMM_WORLD,ierr)
             if(ierr /= 0) then 
                call cq_abort(' ERROR!  MPI_ssend in collect_result ' ,nnd_rem,iprim)
             endif
             !write(io_lun,*) 'MPI_ssend in collect_result END ' &
             !           ,mynode,nnd_rem,iprim
          endif
          !deallocate(recv_array,STAT=stat)
          !if(stat/=0) call cq_abort("Error deallocating recv_array in collect_result: ",nsize)
       enddo ! loop over sending nodes
    endif ! if there are sending nodes

    ! write(io_lun,*) ' MPI_WAIT start for Node',mynode
    !Check whether we have finished all MPI_irecv
    if(comBG%no_recv_node(iprim) > 0) then
       do inode=1,comBG%no_recv_node(iprim)
          nnd_rem=comBG%list_recv_node(inode,iprim)
          if(nnd_rem /= mynode) then
             call MPI_WAIT(nrecv_req(inode),nwait_stat,ierr)
             if(ierr /= 0) then
                call cq_abort('ERROR in MPI_WAIT in collect_result',ierr)
             endif
          endif      ! if it is remote...
       enddo      ! Loop over remote nodes
       call start_timer(tmr_std_allocation)
       deallocate(nrecv_req, STAT=ierr)
       if(ierr/=0) call cq_abort("Error deallocating nrecv_req in collect_result: ",comBG%no_recv_node(iprim),ierr)
       call reg_dealloc_mem(area_basis,comBG%no_recv_node(iprim),type_int)
       call stop_timer(tmr_std_allocation)
    endif       ! if there are sending nodes for iprim...
    if(comBG%no_send_node(iprim)>0) then
       if(allocated(recv_array)) then
          call start_timer(tmr_std_allocation)
          nullify(recv_ptr)
          deallocate(recv_array,STAT=stat)
          if(stat/=0) call cq_abort("Error deallocating recv_array in collect_result: ",nsize)
          call reg_dealloc_mem(area_basis,msize,type_dbl)
          call stop_timer(tmr_std_allocation)
       else
          if(myid==0) write(io_lun,fmt='(2x,"Possible problem in collect_result: recv_array not allocated")')
       end if
    end if
    call my_barrier() ! this is not needed, I think.
    return
  end subroutine collect_result

!!****f* blip_grid_transform_module/inverse_blip_to_grad_new *
!!
!!  NAME 
!!   inverse_blip_to_grad_new
!!  USAGE
!! 
!!  PURPOSE
!!   Takes a gradient wrt support functions on the grid and 
!!   returns gradient wrt blip coefficients
!!  INPUTS
!! 
!! 
!!  USES
!!   do_inverse_blip_to_grad_new
!!   collect_result
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!!   16/01/2001
!!  MODIFICATION HISTORY
!!   17/05/2001 dave
!!    Added ROBODoc, removed subroutine calls
!!   11/06/2001 dave
!!    Added GenComms
!!   2011/11/15 17:00 dave
!!    Removed blip module and variables passed to transform
!!  SOURCE
!!
  subroutine inverse_blip_to_grad_new(myid, direction, dsupport, data_dblip, n_prim)

    use datatypes
    use blip, ONLY: blip_info
    use numbers
    use primary_module, ONLY:bundle
    use set_blipgrid_module, ONLY:naba_blk_supp,comBG
    use block_module,        ONLY: n_pts_in_block  ! = blocks%nm_group(:)
    use GenComms, ONLY: my_barrier
    use species_module, ONLY: nsf_species
    use support_spec_format, ONLY: support_function

    implicit none

    integer,intent(in) :: myid,direction,dsupport,n_prim
    type(support_function),intent(out):: data_dblip(n_prim)

    ! local variables
    integer :: iprim, nsf_recv, spec

    do iprim = 1, bundle%mx_iprim

       !  dsupport shows values on integration grids
       !  Bundle responsible nodes accumulates their 
       ! responsible values on integration grids from
       ! Domain responsible nodes and put them into 
       ! send_array. ( this array is defined in comm_array_module.)

       if(iprim<=bundle%n_prim) then 
          spec = bundle%species(iprim)
          nsf_recv = nsf_species(spec)
       else
          nsf_recv = 0
       end if
       call collect_result (myid,iprim,nsf_recv,dsupport)

       call my_barrier()
       if(iprim <= bundle%n_prim) then
          call do_inverse_blip_to_grad_new (myid, direction, iprim, data_dblip(iprim), nsf_recv)
       endif
    end do
    return
  end subroutine inverse_blip_to_grad_new
!!***

!!****f* blip_grad_transform_module/do_inverse_blip_to_grad_new *
!!
!!  NAME 
!!   do_inverse_blip_to_grad_new
!!  USAGE
!! 
!!  PURPOSE
!!   Does the transform
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!!   23/06/2000
!!  MODIFICATION HISTORY
!!   17/05/2001 dave
!!    Removed the NSF=4 dependencies
!!   2008/06/10 ast
!!    Added timers
!!  SOURCE
!!
  subroutine do_inverse_blip_to_grad_new( myid, direction, iprim,data_dblip, NSF)

    use datatypes
    use numbers
    use global_module,       ONLY:rcellx,rcelly,rcellz, area_basis
    use group_module,        ONLY:blocks
    use primary_module,      ONLY:bundle
    use cover_module,        ONLY:BCS_blocks
    use set_blipgrid_module, ONLY:comBG,naba_blk_supp
    use comm_array_module,   ONLY:send_array

    use block_module,ONLY: n_pts_in_block,nx_in_block,ny_in_block,nz_in_block
    use GenComms, ONLY: cq_abort
    use support_spec_format, ONLY: support_function
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use blip, ONLY: blip_info

    ! Shared variable
    implicit none

    integer,intent(in):: iprim,myid,direction
    integer,intent(in):: nsf
    type(support_function),intent(out):: data_dblip

    ! Local variables
    real(double),allocatable:: inter_1(:,:,:,:)
    real(double),allocatable:: inter_2(:,:,:,:)
    real(double):: rec_sgs, dxxx, dyyy, dzzz
    integer     :: x, y, z, bx, by, bz, imin, imax, nsf1, spec
    integer     :: nxmin_grid,nxmax_grid
    integer     :: nymin_grid,nymax_grid
    integer     :: nzmin_grid,nzmax_grid
    integer, allocatable, dimension(:) :: imin_for_z, imax_for_z

    real(double), allocatable, dimension(:) :: splines
    real(double), allocatable, dimension(:,:) :: splines_for_z
    real(double):: ys, dsum(NSF), yss

    integer     :: ierr,stat
    integer     :: ix_grid,iy_grid,iz_grid,ix,iy,iz
    real(double):: dcellx_block,dcelly_block,dcellz_block
    real(double):: dcellx_grid, dcelly_grid, dcellz_grid
    integer     :: ind_blip,igrid,ind_blk,nx_blk,ny_blk,nz_blk
    integer     :: ncover_yz,naba_blk

    !!   nblkx,nblky,nblkz : number of integration grid points in a block
    !!     this variables should be passed from (blocks)
    integer :: nblkx,nblky,nblkz

    spec = bundle%species(iprim)
    call start_timer(tmr_std_allocation)
    allocate(inter_1(NSF,-blip_info(spec)%BlipArraySize:blip_info(spec)%BlipArraySize,&
         -blip_info(spec)%BlipArraySize:blip_info(spec)%BlipArraySize,2*blip_info(spec)%Extent+1),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to inter_1 !',stat)
    call reg_alloc_mem(area_basis,nsf*(2*blip_info(spec)%BlipArraySize+1)*(2*blip_info(spec)%BlipArraySize+1)*&
         (2*blip_info(spec)%Extent+1),type_dbl)
    call stop_timer(tmr_std_allocation)
    inter_1=zero
    call start_timer(tmr_std_allocation)
    allocate(inter_2(NSF,-blip_info(spec)%BlipArraySize:blip_info(spec)%BlipArraySize,&
         2*blip_info(spec)%Extent+1,2*blip_info(spec)%Extent+1),STAT=stat)
    if(stat/=0) call cq_abort('Error allocating memory to inter_2 !',stat)
    call reg_alloc_mem(area_basis,nsf*(2*blip_info(spec)%BlipArraySize+1)*&
         (2*blip_info(spec)%Extent+1)*(2*blip_info(spec)%Extent+1),type_dbl)
    allocate(splines( -blip_info(spec)%BlipArraySize:blip_info(spec)%BlipArraySize ), &
         splines_for_z(-blip_info(spec)%BlipArraySize:blip_info(spec)%BlipArraySize,1:2*blip_info(spec)%Extent+1))
    call reg_alloc_mem(area_basis,(2*blip_info(spec)%BlipArraySize+1)*(2*blip_info(spec)%Extent+2),type_dbl)
    allocate(imin_for_z(2*blip_info(spec)%Extent+1),imax_for_z(2*blip_info(spec)%Extent+1),STAT=stat)
    call reg_alloc_mem(area_basis,2*(2*blip_info(spec)%Extent+1),type_dbl)
    call stop_timer(tmr_std_allocation)
    inter_2=zero
    imin_for_z = 0
    imax_for_z = 0
    splines = zero
    splines_for_z = zero

    rec_sgs = one / blip_info(spec)%SupportGridSpacing
    nblkx=nx_in_block
    nblky=ny_in_block
    nblkz=nz_in_block

    dcellx_block=rcellx/blocks%ngcellx; dcellx_grid=dcellx_block/nblkx
    dcelly_block=rcelly/blocks%ngcelly; dcelly_grid=dcelly_block/nblky
    dcellz_block=rcellz/blocks%ngcellz; dcellz_grid=dcellz_block/nblkz

    nxmin_grid=(naba_blk_supp%nxmin(iprim)-1)*nblkx
    nxmax_grid= naba_blk_supp%nxmax(iprim)*nblkx-1
    nymin_grid=(naba_blk_supp%nymin(iprim)-1)*nblky
    nymax_grid= naba_blk_supp%nymax(iprim)*nblky-1
    nzmin_grid=(naba_blk_supp%nzmin(iprim)-1)*nblkz
    nzmax_grid= naba_blk_supp%nzmax(iprim)*nblkz-1

    !check extent
    ierr=0
    if(2*blip_info(spec)%Extent < nxmax_grid-nxmin_grid) then
       write(io_lun,*) myid,' ERROR in do_inverse_blip : extent for x dir ::', &
            ' extent nxmin_grid nxmax_grid = ',blip_info(spec)%Extent,nxmin_grid,nxmax_grid
       ierr=1
    endif
    if(2*blip_info(spec)%Extent < nymax_grid-nymin_grid) then
       write(io_lun,*) myid,' ERROR in do_inverse_blip : extent for y dir ::', &
            ' extent nymin_grid nymax_grid = ',blip_info(spec)%Extent,nymin_grid,nymax_grid
       ierr=1
    endif
    if(2*blip_info(spec)%Extent < nzmax_grid-nzmin_grid) then
       write(io_lun,*) myid,' ERROR in do_inverse_blip : extent for z dir ::', &
            ' extent nymin_grid nymax_grid = ',blip_info(spec)%Extent,nzmin_grid,nzmax_grid
       ierr=1
    endif
    if(ierr /= 0) call cq_abort('do_inverse_blip ',ierr)
    DO z = nzmin_grid , nzmax_grid
       iz_grid=z-nzmin_grid+1
       dzzz=z*dcellz_grid-bundle%zprim(iprim)
       imin_for_z(iz_grid)=MAX(-blip_info(spec)%BlipArraySize, &
            int((dzzz-blip_info(spec)%BlipWidth*half)*rec_sgs))
       imax_for_z(iz_grid)=MIN( blip_info(spec)%BlipArraySize, &
            int((dzzz+blip_info(spec)%BlipWidth*half)*rec_sgs))
       if(imin_for_z(iz_grid) > blip_info(spec)%BlipArraySize) then
          imin_for_z(iz_grid)=blip_info(spec)%BlipArraySize
       endif
       if(imax_for_z(iz_grid) < -blip_info(spec)%BlipArraySize) then
          imax_for_z(iz_grid)=-blip_info(spec)%BlipArraySize
       endif

       DO bz = imin_for_z(iz_grid), imax_for_z(iz_grid)
          yss = dzzz - bz * blip_info(spec)%SupportGridSpacing 
          ys = ABS( yss )
          ys = ys * blip_info(spec)%FourOnBlipWidth

          if (yss.ge.zero) then
             yss = blip_info(spec)%FourOnBlipWidth
          else
             yss = -blip_info(spec)%FourOnBlipWidth
          end if

          IF( ys .LE. one ) THEN
             if (direction.ne.3) then
                splines_for_z(bz,iz_grid) = one +                      &
                     ys * ys * ( 0.75_double * ys - 1.5_double )
             else
                splines_for_z(bz,iz_grid) = yss * ys * (2.25_double * ys - 3.0_double )
             end if

          ELSE IF( ys .LE. two ) THEN
             if (direction.ne.3) then
                splines_for_z(bz,iz_grid) = 0.25_double * ( two - ys ) *    &
                     ( two - ys ) * ( two - ys )
             else
                splines_for_z(bz,iz_grid) = &
                     -yss * 0.75_double * ( two - ys ) * ( two - ys )
             end if
          ELSE
             splines_for_z(bz,iz_grid) = zero
          END IF
       END DO
    ENDDO

    ncover_yz=BCS_blocks%ncovery*BCS_blocks%ncoverz
    igrid=0
    DO naba_blk=1,naba_blk_supp%no_naba_blk(iprim)! naba blks (in NOPG order)
       ! iprim : prim. seq. # of atm
       ind_blk=naba_blk_supp%list_naba_blk(naba_blk,iprim) !CC in BCS_blocks
       ind_blk=ind_blk-1

       nx_blk=ind_blk/ncover_yz
       ny_blk=(ind_blk-nx_blk*ncover_yz)/BCS_blocks%ncoverz
       nz_blk=ind_blk-nx_blk*ncover_yz-ny_blk*BCS_blocks%ncoverz

       nx_blk=nx_blk-BCS_blocks%nspanlx+BCS_blocks%nx_origin
       ny_blk=ny_blk-BCS_blocks%nspanly+BCS_blocks%ny_origin
       nz_blk=nz_blk-BCS_blocks%nspanlz+BCS_blocks%nz_origin

       do iz=1,nblkz    ! z-direction in a block
          z=nblkz*(nz_blk-1)+iz-1
          iz_grid=z-nzmin_grid+1

          imin=imin_for_z(iz_grid)
          imax=imax_for_z(iz_grid)
          do iy=1,nblky    ! y-direction in a block
             y=nblky*(ny_blk-1)+iy-1
             iy_grid=y-nymin_grid+1
             if(iy_grid < 1 .OR. iy_grid > nymax_grid-nymin_grid+1) then
                call cq_abort(' ERROR iy_grid = ',iy_grid,nymax_grid-nymin_grid+1)
             endif
             do ix=1,nblkx    ! x-direction in a block
                x=nblkx*(nx_blk-1)+ix-1
                ix_grid=x-nxmin_grid+1
                if(ix_grid < 1 .OR. ix_grid > nxmax_grid-nxmin_grid+1) then
                   call cq_abort(' ERROR ix_grid = ',ix_grid,nxmax_grid-nxmin_grid+1)
                endif

                igrid=igrid+1  ! seq. no. of integration grids

                if(igrid > naba_blk_supp%no_naba_blk(iprim)*n_pts_in_block .OR. igrid < 1) then
                   call cq_abort('ERROR!!!!!!! ',igrid,naba_blk_supp%no_naba_blk(iprim)*n_pts_in_block)
                endif
                do nsf1 = 1,NSF ! data_naba_grid(1,igrid)
                   dsum(nsf1) = send_array(NSF*igrid-NSF+nsf1) 
                enddo
                DO bz = imin,imax
                   ys = splines_for_z(bz,iz_grid)
                   do nsf1 = 1,NSF
                      inter_2(nsf1,bz,ix_grid,iy_grid)= &
                           inter_2( nsf1,bz,ix_grid,iy_grid)+dsum(nsf1)*ys
                   enddo
                END DO
             enddo    ! x-direction in a block
          enddo    ! y-direction in a block
       enddo    ! z-direction in a block
    END DO    ! Loop over naba blocks for this primary atom
    DO y = nymin_grid, nymax_grid
       iy_grid=y-nymin_grid+1
       dyyy=y*dcelly_grid-bundle%yprim(iprim)
       imin = MAX(-blip_info(spec)%BlipArraySize,int((dyyy-blip_info(spec)%BlipWidth*half)*rec_sgs))
       imax = MIN( blip_info(spec)%BlipArraySize,int((dyyy+blip_info(spec)%BlipWidth*half)*rec_sgs))
       if(imin >  blip_info(spec)%BlipArraySize) imin= blip_info(spec)%BlipArraySize
       if(imax < -blip_info(spec)%BlipArraySize) imax=-blip_info(spec)%BlipArraySize
       DO by = imin, imax
          yss = dyyy - by * blip_info(spec)%SupportGridSpacing 
          ys = ABS( yss )
          ys = ys * blip_info(spec)%FourOnBlipWidth

          if (yss.ge.zero) then
             yss = blip_info(spec)%FourOnBlipWidth
          else
             yss =-blip_info(spec)%FourOnBlipWidth
          end if

          IF( ys .LE. one ) THEN
             if (direction.ne.2) then
                splines( by ) = one +  &
                     ys * ys * ( 0.75_double * ys - 1.5_double )
             else
                splines( by ) = yss * ys * (2.25_double * ys - 3.0_double )
             end if
          ELSE IF( ys .LE. two ) THEN
             if (direction.ne.2) then
                splines( by ) = 0.25_double * ( two - ys ) *  &
                     ( two - ys ) * ( two - ys )
             else
                splines( by ) = -yss * 0.75_double * ( two - ys ) * ( two - ys )
             end if
          ELSE
             splines( by ) = zero
          END IF
       END DO
       DO x = nxmin_grid, nxmax_grid
          ix_grid=x-nxmin_grid+1
          DO bz = -blip_info(spec)%BlipArraySize, blip_info(spec)%BlipArraySize
             do nsf1 = 1,NSF
                dsum(nsf1) = inter_2(nsf1,bz,ix_grid,iy_grid )
             enddo
             DO by = imin, imax
                ys = splines( by )
                do nsf1 = 1,NSF
                   inter_1(nsf1,by,bz,ix_grid)=dsum(nsf1)*ys + &
                        inter_1(nsf1,by,bz,ix_grid)
                enddo
             END DO
          END DO
       END DO
    END DO
    DO x = nxmin_grid, nxmax_grid
       ix_grid=x-nxmin_grid+1
       dxxx=x*dcellx_grid-bundle%xprim(iprim)
       imin = MAX(-blip_info(spec)%BlipArraySize,int((dxxx-blip_info(spec)%BlipWidth*half)*rec_sgs))
       imax = MIN( blip_info(spec)%BlipArraySize,int((dxxx+blip_info(spec)%BlipWidth*half)*rec_sgs))
       if(imin >  blip_info(spec)%BlipArraySize) imin= blip_info(spec)%BlipArraySize
       if(imax < -blip_info(spec)%BlipArraySize) imax=-blip_info(spec)%BlipArraySize
       DO bx =imin, imax
          yss = dxxx - bx * blip_info(spec)%SupportGridSpacing
          ys = ABS( yss )
          ys = ys * blip_info(spec)%FourOnBlipWidth

          if (yss.ge.zero) then
             yss = blip_info(spec)%FourOnBlipWidth
          else
             yss =-blip_info(spec)%FourOnBlipWidth
          end if

          IF( ys .LE. one ) THEN
             if (direction.ne.1) then
                splines( bx ) = one + &
                     ys * ys * ( 0.75_double * ys - 1.5_double )
             else
                splines( bx ) = yss * ys * (2.25_double * ys - 3.0_double )
             end if
          ELSE IF( ys .LE. two ) THEN
             if (direction.ne.1) then
                splines( bx ) = 0.25_double * ( two - ys ) *  &
                     ( two - ys ) * ( two - ys )
             else
                splines( bx ) = -yss * 0.75_double * ( two - ys ) * ( two - ys )
             end if
          ELSE
             splines( bx ) = zero
          END IF
       END DO
       DO bz = -blip_info(spec)%BlipArraySize, blip_info(spec)%BlipArraySize
          DO by = -blip_info(spec)%BlipArraySize, blip_info(spec)%BlipArraySize
             do nsf1 = 1,NSF
                dsum(nsf1) = inter_1( nsf1,by,bz,ix_grid)
             enddo
             DO bx = imin, imax
                ind_blip = blip_info(spec)%blip_number(bx,by,bz)
                IF(ind_blip /= 0 ) THEN
                   ys = splines( bx )
                   do nsf1 = 1,NSF
                      data_dblip%supp_func(nsf1)%coefficients(ind_blip)=dsum(nsf1)*ys+ &
                           data_dblip%supp_func(nsf1)%coefficients(ind_blip)
                   enddo
                END IF
             END DO
          END DO
       END DO
    END DO
    call start_timer(tmr_std_allocation)
    deallocate(inter_1,inter_2,STAT=stat)
    if(stat/=0) call cq_abort('Error deallocating memory to inter_1 !',stat)
    call reg_dealloc_mem(area_basis,nsf*(2*blip_info(spec)%BlipArraySize+1)*&
         (2*blip_info(spec)%BlipArraySize+1)*(2*blip_info(spec)%Extent+1),type_dbl)
    call reg_dealloc_mem(area_basis,nsf*(2*blip_info(spec)%BlipArraySize+1)*&
         (2*blip_info(spec)%Extent+1)*(2*blip_info(spec)%Extent+1),type_dbl)
    call reg_dealloc_mem(area_basis,size(send_array),type_dbl)
    deallocate(splines,splines_for_z,STAT=stat)
    call reg_dealloc_mem(area_basis,(2*blip_info(spec)%BlipArraySize+1)*(2*blip_info(spec)%Extent+2),type_dbl)
    if(stat/=0) call cq_abort('Error deallocating memory to splines !',stat)

    deallocate(send_array,STAT=stat)
    if(stat/=0) call cq_abort('Error deallocating send_array in inverse_blip_to_grad !',stat)
    deallocate(imin_for_z,imax_for_z,STAT=stat)
    call reg_dealloc_mem(area_basis,2*(2*blip_info(spec)%Extent+1),type_dbl)
    call stop_timer(tmr_std_allocation)
    return
  end subroutine do_inverse_blip_to_grad_new
!!***
end module blip_grid_transform_module
