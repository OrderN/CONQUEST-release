! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module set_blipgrid_module
! ------------------------------------------------------------------------------
! Code area 8: grids and indexing
! ------------------------------------------------------------------------------

!!****h* Conquest/set_blipgrid_module
!!  NAME
!!   set_blipgrid_module
!!  PURPOSE
!!   Sets up the relationships between the integration grid and the atoms - in 
!!   particular (and historically) between the blips on atoms and the integration grid
!!  USES
!!
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!!   Mid-2000
!!  MODIFICATION HISTORY
!!   12:03, 25/09/2002 mjg & drb 
!!    Added headers, and created neighbour lists for atomic densities
!!   08:19, 2003/06/11 dave
!!    Added ROBOdoc headers to about half routines and TM's debug and tidy statements
!!   2006/07/05 08:14 dave
!!    Various changes for variable NSF, and consolidation of different local orbital definitions
!!   2008/02/06 08:35 dave
!!    Changed for output to file not stdout
!!   2008/05/16 ast
!!    Added timers
!!   2014/09/15 18:30 lat
!!    fixed call start/stop_timer to timer_module (not timer_stdlocks_module !)
!!   2016/07/06 17:30 nakata
!!    Renamed comm_in_BG -> comm_in_BtoG, comBG -> comm_naba_blocks_of_atoms and mx_func_BG -> mx_func_BtoG
!!   2016/07/14 16:30 nakata
!!    Renamed naba_blk_supp -> naba_blocks_of_atoms
!!   2016/07/20 16:30 nakata
!!    Renamed naba_atm -> naba_atoms_of_blocks, halo_atm -> halo_atoms_of_blocks
!!   2016/08/01 17:30 nakata
!!    Introduced atomf
!!  SOURCE
!!
module set_blipgrid_module

  use global_module,   only: sf, nlpf, paof, dens, atomf, iprint_index, io_lun
  use naba_blk_module, only: &
       comm_in_BtoG, naba_blk_of_atm, naba_atm_of_blk, halo_atm_of_blk, &
       alloc_comm_in_BtoG1, alloc_comm_in_BtoG2, alloc_naba_blk, alloc_naba_atm, alloc_halo_atm

  use comm_array_module,      only: isend_array, irecv_array
  use timer_module,           only: start_timer, stop_timer
  use timer_stdclocks_module, only: tmr_std_indexing, tmr_std_allocation

  !check_commarray_int!,check_commarray_real

  implicit none

  ! Definition 
  integer, parameter    :: mx_func_BtoG=4
  !integer, parameter    :: supp=1
  !integer, parameter    :: proj=2
  !integer, parameter    :: dens=3

  type(comm_in_BtoG)           :: comm_naba_blocks_of_atoms
  type(naba_blk_of_atm)        :: naba_blocks_of_atoms
  type(naba_atm_of_blk),target :: naba_atoms_of_blocks(mx_func_BtoG)
  type(halo_atm_of_blk),target :: halo_atoms_of_blocks(mx_func_BtoG)

  ! pointer used in recv_array_BtoG and make_table_BtoG
  integer,pointer :: naba_blk_rem(:)
  integer,pointer :: offset_rem(:)

  ! RCS tag for object file identification
  character(len=80), save, private :: RCSid = "$Id$"
!!***

contains

!!****f* set_blipgrid_module/set_blipgrid *
!!
!!  NAME 
!!   set_blipgrid
!!  USAGE
!! 
!!  PURPOSE
!!    main subroutine in this module
!!    set_blipgrid
!!     |--  alloc_comm_in_BtoG
!!     |-- (alloc_naba_blk,alloc_naba_atm,alloc_halo_atm)
!!     |                                  in naba_blk_module
!!     |--  get_naba_BCSblk
!!     |---  dist_blk_atom
!!     |--  get_naba_DCSprt
!!     |--  make_sendinfo_BtoG
!!     |--  make_table_BtoG
!!              |-- send_array_BtoG
!!              |-- recv_array_BtoG
!!              |-- make_table
!! 
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!!   31/5/2000
!!  MODIFICATION HISTORY
!!   8/9/2000 TM
!!    I added subroutines in make_table_BtoG.f90
!!   12:16, 25/09/2002 mjg & drb 
!!    Added use atomic_density to get rcut_dens
!!   2006/07/06 08:27 dave
!!    Various changes including use of cq_abort
!!   2008/05/16 ast
!!    Added timers
!!   2011/11/17 10:41 dave
!!    Changes for new blip data (principally only defining Extent in this routine)
!!   2016/07/06 17:30 nakata
!!    Renamed comm_in_BG -> comm_in_BtoG, comBG -> comm_naba_blocks_of_atoms and max_recv_node_BG -> max_recv_node_BtoG
!!   2016/07/14 16:30 nakata
!!    Renamed naba_blk_supp -> naba_blocks_of_atoms
!!   2016/07/20 16:30 nakata
!!    Renamed naba_atm -> naba_atoms_of_blocks, halo_atm -> halo_atoms_of_blocks
!!   2016/08/01 17:30 nakata
!!    Introduced atomf instead of sf and paof
!!   2016/08/08 15:30 nakata
!!    Renamed rcut_supp -> rcut_atomf
!!  SOURCE
!!
  subroutine set_blipgrid(myid,rcut_atomf,rcut_proj)
    use datatypes
    use global_module,ONLY:numprocs, sf, nlpf, paof, atomf
    use primary_module,ONLY: bundle, domain
    use cover_module,  ONLY: DCS_parts
    use atomic_density, ONLY: rcut_dens
    use GenComms, ONLY: cq_abort
    use functions_on_grid, ONLY: gridsize
    use block_module, ONLY: n_pts_in_block, nx_in_block,ny_in_block,nz_in_block
    use blip, ONLY: blip_info
    use species_module, ONLY: n_species
    !use pseudopotential_common
    use dimens, ONLY: n_grid_x, n_grid_y, n_grid_z, r_super_x, r_super_y, r_super_z
    
    implicit none
    integer,intent(in)      ::myid
    real(double),intent(in), dimension(n_species) ::rcut_atomf,rcut_proj

    integer :: iprim,max_naba_blk,max_naba_atm,max_naba_part, max_halo_part
    integer :: nxmin_grid, nxmax_grid,  nymin_grid, nymax_grid, nzmin_grid, nzmax_grid, xextent, yextent, zextent
    integer :: max_recv_node_BtoG, max_send_node, max_sent_pairs, max_recv_call
    integer :: thisextent,i, spec
    logical :: warn

    call start_timer(tmr_std_indexing)
    ! Find maxima and allocate derived types
    call get_naba_BCSblk_max(rcut_atomf, max_naba_blk, max_recv_node_BtoG)
    call alloc_comm_in_BtoG1(comm_naba_blocks_of_atoms,bundle%mx_iprim,max_recv_node_BtoG)
    !call alloc_comm_in_BtoG(comm_naba_blocks_of_atoms,bundle%mx_iprim,mx_recv_node_BtoG,mx_send_node_BtoG,&
    !     mx_sent_pair_BtoG,mx_recv_call_BtoG)
    call alloc_naba_blk(naba_blocks_of_atoms,bundle%mx_iprim,max_naba_blk)
    ! Allocate naba_atoms_of_blocks derived types
!!! nakata2
    call get_naba_DCSprt_max(rcut_atomf,max_naba_part,max_naba_atm,max_halo_part, max_recv_node_BtoG)
    call alloc_naba_atm(naba_atoms_of_blocks(atomf), domain%mx_ngonn,max_naba_part,max_naba_atm)
    naba_atoms_of_blocks(atomf)%function_type = atomf
    call alloc_halo_atm(halo_atoms_of_blocks(atomf), max_recv_node_BtoG,max_halo_part,DCS_parts%mx_mcover)
!!! nakata2 --- delete from here later
!    ! Support and PAOs have same cutoff and maxima
!    call get_naba_DCSprt_max(rcut_supp,max_naba_part,max_naba_atm,max_halo_part, max_recv_node_BtoG)
!    call alloc_naba_atm(naba_atoms_of_blocks(sf), domain%mx_ngonn,max_naba_part,max_naba_atm)
!    naba_atoms_of_blocks(sf)%function_type = sf
!    call alloc_halo_atm(halo_atoms_of_blocks(sf), max_recv_node_BtoG,max_halo_part,DCS_parts%mx_mcover)
!    ! PAOs
!    call alloc_naba_atm(naba_atoms_of_blocks(paof), domain%mx_ngonn,max_naba_part,max_naba_atm)
!    naba_atoms_of_blocks(paof)%function_type = paof
!    call alloc_halo_atm(halo_atoms_of_blocks(paof), max_recv_node_BtoG,max_halo_part,DCS_parts%mx_mcover)
!!! nakata2 --- delet up to here later
!!! nakata2 end
    ! Projectors
    call get_naba_DCSprt_max(rcut_proj,max_naba_part,max_naba_atm,max_halo_part, max_recv_node_BtoG)
    call alloc_naba_atm(naba_atoms_of_blocks(nlpf), domain%mx_ngonn,max_naba_part,max_naba_atm)
    naba_atoms_of_blocks(nlpf)%function_type = nlpf
    warn = .false.
    do i=1,n_species
       if(rcut_atomf(i) < rcut_proj(i)) warn = .true.
    end do
    if(.NOT.warn) then
       call alloc_halo_atm(halo_atoms_of_blocks(nlpf), max_recv_node_BtoG,max_halo_part,DCS_parts%mx_mcover)
    else
       ! Note: this is a problem because max_recv_node_BtoG isn't defined
       if(myid==0) write(io_lun,*) ' WARNING !! '
       if(myid==0) write(io_lun,*) ' before alloc_halo_atm : rcut_atomf < rcut_proj '
       call alloc_halo_atm(halo_atoms_of_blocks(nlpf), numprocs,max_halo_part,DCS_parts%mx_mcover)
    endif
    ! Densities
    call get_naba_DCSprt_max(rcut_dens,max_naba_part,max_naba_atm,max_halo_part, max_recv_node_BtoG)
    call alloc_naba_atm(naba_atoms_of_blocks(dens), domain%mx_ngonn,max_naba_part,max_naba_atm)
    naba_atoms_of_blocks(dens)%function_type = dens
    call alloc_halo_atm(halo_atoms_of_blocks(dens), max_recv_node_BtoG,max_halo_part,DCS_parts%mx_mcover)


    !make lists of neighbour blocks of primary atoms and
    ! information for sending blip-grid transformed support functions
    call get_naba_BCSblk(rcut_atomf,naba_blocks_of_atoms,comm_naba_blocks_of_atoms)
    ! Calculate likely extents
    do spec = 1,n_species
       xextent = int((rcut_atomf(spec)*n_grid_x/r_super_x)+0.5)
       yextent = int((rcut_atomf(spec)*n_grid_y/r_super_y)+0.5)
       zextent = int((rcut_atomf(spec)*n_grid_z/r_super_z)+0.5)
       blip_info(spec)%Extent = MAX(xextent,MAX(yextent,zextent))
    end do
    do iprim=1,bundle%n_prim
    !do iprim=1,bundle%mx_iprim
       thisextent = 0
       if(iprim == 1) then
          max_naba_blk=naba_blocks_of_atoms%no_naba_blk(iprim)
       elseif ( naba_blocks_of_atoms%no_naba_blk(iprim) > max_naba_blk) then
          max_naba_blk=naba_blocks_of_atoms%no_naba_blk(iprim)
       endif
       ! Calculate extent
       nxmin_grid=(naba_blocks_of_atoms%nxmin(iprim)-1)*nx_in_block  
       nxmax_grid= naba_blocks_of_atoms%nxmax(iprim)*nx_in_block-1
       nymin_grid=(naba_blocks_of_atoms%nymin(iprim)-1)*ny_in_block  
       nymax_grid= naba_blocks_of_atoms%nymax(iprim)*ny_in_block-1
       nzmin_grid=(naba_blocks_of_atoms%nzmin(iprim)-1)*nz_in_block  
       nzmax_grid= naba_blocks_of_atoms%nzmax(iprim)*nz_in_block-1
       if(iprim==1) then
          xextent=(nxmax_grid-nxmin_grid+1)/2
          yextent=(nymax_grid-nymin_grid+1)/2
          zextent=(nzmax_grid-nzmin_grid+1)/2
       else
          xextent=max(xextent,(nxmax_grid-nxmin_grid+1)/2)
          yextent=max(yextent,(nymax_grid-nymin_grid+1)/2)
          zextent=max(zextent,(nzmax_grid-nzmin_grid+1)/2)
       end if
       thisextent = max(zextent,max(yextent,xextent))
       spec = bundle%species(iprim)
       blip_info(spec)%Extent = max(thisextent,blip_info(spec)%Extent)
    enddo
    !make lists of neighbour and halo atoms of primary blocks
    ! for support and projector functions
!!! nakata2
    call get_naba_DCSprt(rcut_atomf,naba_atoms_of_blocks(atomf),halo_atoms_of_blocks(atomf))
    !do iprim_blk=1,domain%groups_on_node
    gridsize(atomf) = 0
    do iprim=1,domain%mx_ngonn
       gridsize(atomf) = gridsize(atomf) + naba_atoms_of_blocks(atomf)%no_of_orb(iprim)*n_pts_in_block
    end do
!!! nakata2 --- delete from here later
!    call get_naba_DCSprt(rcut_supp,naba_atoms_of_blocks(sf),halo_atoms_of_blocks(sf))
!    !do iprim_blk=1,domain%groups_on_node
!    gridsize(sf) = 0
!    do iprim=1,domain%mx_ngonn
!       gridsize(sf) = gridsize(sf) + naba_atoms_of_blocks(sf)%no_of_orb(iprim)*n_pts_in_block
!    end do
!    call get_naba_DCSprt(rcut_supp,naba_atoms_of_blocks(paof),halo_atoms_of_blocks(paof))
!    gridsize(paof) = 0
!    do iprim=1,domain%mx_ngonn
!       gridsize(paof) = gridsize(paof) + naba_atoms_of_blocks(paof)%no_of_orb(iprim)*n_pts_in_block
!    end do
!!! nakata2 --- delete up to here later
!!! nakata2 end
    call get_naba_DCSprt(rcut_proj,naba_atoms_of_blocks(nlpf),halo_atoms_of_blocks(nlpf))
    !do iprim_blk=1,domain%groups_on_node
    gridsize(nlpf) = 0
    do iprim=1,domain%mx_ngonn
       gridsize(nlpf) = gridsize(nlpf) + naba_atoms_of_blocks(nlpf)%no_of_orb(iprim)*n_pts_in_block
    end do
    !write(io_lun,*) 'gridsizes (atomf, nlpf): ',gridsize(atomf),gridsize(nlpf)
    ! Now add atomic density tables
    call get_naba_DCSprt(rcut_dens,naba_atoms_of_blocks(dens),halo_atoms_of_blocks(dens))
    ! End of atomic density tables

    !prepare for receiving
    call make_sendinfo_BtoG_max(myid,naba_atoms_of_blocks(atomf),max_send_node,max_sent_pairs,max_recv_call)
    call alloc_comm_in_BtoG2(comm_naba_blocks_of_atoms,max_send_node, max_sent_pairs,max_recv_call)
    call make_sendinfo_BtoG(myid,naba_atoms_of_blocks(atomf),comm_naba_blocks_of_atoms)
    !making tables showing where to put sent BtoG transformed atomic func. (SF or PAOs)
    ! Now make_table_BtoG is in <make_table_BtoG.f90>.
    call make_table_BtoG(myid,comm_naba_blocks_of_atoms)
    call stop_timer(tmr_std_indexing)
    return
  end subroutine set_blipgrid
!!***

!!****f* set_blipgrid_module/get_naba_BCSblk *
!!
!!  NAME 
!!   get_naba_BCSblk
!!  USAGE
!! 
!!  PURPOSE
!!   finds neighbour blocks for each of my primary atoms
!!   from my BCS_blk covering sets
!!
!! <STRUCTURE>
!!     do np= primary set (NOPG order)
!!      do ni= atoms in the np-th partition
!!       do iblock= BCS_blocks (Bundle Cover Set blocks)
!!          call distsq_blk_atom
!!         if(neighbour) then
!!           makes naba_blk (type : naba_blk_of_atm) including the pair's offset
!!            and  some parts (information of sending of BtoG transforms) 
!!                  of comm_naba_blocks_of_atoms(type: comm_in_BtoG)
!!         endif
!!       enddo
!!      enddo
!!     enddo
!! 
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!! 
!!  MODIFICATION HISTORY
!!   2006/07/06 08:28 dave
!!    Added cq_abort
!!   2016/07/06 17:30 nakata
!!    Renamed comm_in_BG -> comm_in_BtoG and comBG -> comm_naba_blocks_of_atoms
!!  SOURCE
!!
  subroutine get_naba_BCSblk(rcut,naba_blk,comm_naba_blocks_of_atoms)

    use datatypes
    use numbers,        ONLY: very_small, pi,three,four
    use group_module,   ONLY: blocks,parts
    use primary_module, ONLY: bundle
    use cover_module,   ONLY: BCS_blocks,DCS_parts
    use global_module,  ONLY: rcellx,rcelly,rcellz,numprocs
    use block_module,   ONLY: nx_in_block,ny_in_block,nz_in_block
    use GenComms, ONLY: cq_abort, myid
    use species_module, ONLY: n_species

    implicit none

    !Dummy Arguments
    real(double),intent(in),dimension(n_species) :: rcut
    type(naba_blk_of_atm),intent(inout) ::naba_blk
    type(comm_in_BtoG),intent(inout) ::comm_naba_blocks_of_atoms
    !Local variables
    real(double) :: dcellx_part,dcelly_part,dcellz_part
    real(double),dimension(n_species) :: rcutsq
    real(double) :: dcellx_block,dcelly_block,dcellz_block
    real(double) :: dcellx_grid,dcelly_grid,dcellz_grid ! 4/Aug/2000 TM
    real(double) :: xmin_p,ymin_p,zmin_p
    real(double) :: xmin,xmax,ymin,ymax,zmin,zmax
    real(double) :: xatom,yatom,zatom,distsq
    integer :: inp,ncoverz,ncoveryz,np,nx,ny,nz,ni,iblock,ind_block,nnd_rem
    integer :: ncoverz_rem,ncoveryz_rem,spec
    ! 28/Jul/20000
    integer :: ncoverx_add,ncovery_add,ncoverz_add
    integer :: irc,ierr
    integer :: ierror=0,ind_part,nx1,ny1,nz1   ! ?? temporary
    real(double) :: tmp(bundle%n_prim)
    logical :: flag_new

    do ni = 1,n_species
       rcutsq(ni)=rcut(ni)*rcut(ni)
    end do
    inp=0 
    naba_blk%no_naba_blk=0 
    comm_naba_blocks_of_atoms%no_recv_node=0
    comm_naba_blocks_of_atoms%no_naba_blk=0

    ncoverz=BCS_blocks%ncoverz
    ncoveryz=BCS_blocks%ncovery*BCS_blocks%ncoverz

    dcellx_part=rcellx/parts%ngcellx ; dcellx_block=rcellx/blocks%ngcellx
    dcelly_part=rcelly/parts%ngcelly ; dcelly_block=rcelly/blocks%ngcelly
    dcellz_part=rcellz/parts%ngcellz ; dcellz_block=rcellz/blocks%ngcellz

    dcellx_grid=dcellx_block/nx_in_block
    dcelly_grid=dcelly_block/ny_in_block
    dcellz_grid=dcellz_block/nz_in_block

    !- CHECK BCS_blocks ---
    !  write(io_lun,*) ' ng_cover of BCS_blocks ',&
    !         BCS_blocks%ng_cover
    !  write(io_lun,*) ' ncovers of BCS_blocks ',&
    !         BCS_blocks%ncoverx,BCS_blocks%ncovery,BCS_blocks%ncoverz
    !  write(io_lun,*) ' origin of BCS_blocks ',&
    !         BCS_blocks%nx_origin,BCS_blocks%ny_origin,BCS_blocks%nz_origin
    !  write(io_lun,*) ' nleft of BCS_blocks ',&
    !         BCS_blocks%nspanlx,BCS_blocks%nspanly,BCS_blocks%nspanly

    do np=1,bundle%groups_on_node  ! primary partitions in bundle
       if(bundle%nm_nodgroup(np) > 0) then  ! Are there atoms?
          nx=bundle%idisp_primx(np)+bundle%nx_origin
          ny=bundle%idisp_primy(np)+bundle%ny_origin
          nz=bundle%idisp_primz(np)+bundle%nz_origin
          !DEBUG
          nx1=1+mod(nx-1+parts%ngcellx,parts%ngcellx)
          ny1=1+mod(ny-1+parts%ngcelly,parts%ngcelly)
          nz1=1+mod(nz-1+parts%ngcellz,parts%ngcellz)
          ind_part=(nx1-1)*parts%ngcelly*parts%ngcellz+(ny1-1)*parts%ngcellz+nz1
          if(np /= parts%i_cc2seq(ind_part)) write(io_lun,*) 'BUNDLE ERROR!!',&
               ' nx1,ny1,nz1 = ',nx1,ny1,nz1,'ind',ind_part,'np & i_cc2 = ', &
               np,parts%i_cc2seq(ind_part)

          !l.h.s. corner of the partition
          xmin_p= dcellx_part*(nx-1)
          ymin_p= dcelly_part*(ny-1)
          zmin_p= dcellz_part*(nz-1)

          ierror=0
          do ni=1,bundle%nm_nodgroup(np)  ! number of atoms in the partition
             inp=inp+1                      ! inp is a primary number of the atom

             if(inp > bundle%mx_iprim) call cq_abort(' Error : iprim in get_naba_BCSblk ',inp, bundle%mx_iprim)
             xatom=bundle%xprim(inp)
             yatom=bundle%yprim(inp)
             zatom=bundle%zprim(inp)
             spec = bundle%species(inp)
             do iblock=1,BCS_blocks%ng_cover             !NOPG order (cover blks)
                ind_block=BCS_blocks%lab_cover(iblock)-1 !(CC in a cov. set)-1
                nx=ind_block/ncoveryz
                ny=(ind_block-nx*ncoveryz)/ncoverz
                nz=ind_block-nx*ncoveryz-ny*ncoverz

                nx=nx-BCS_blocks%nspanlx+BCS_blocks%nx_origin
                ny=ny-BCS_blocks%nspanly+BCS_blocks%ny_origin
                nz=nz-BCS_blocks%nspanlz+BCS_blocks%nz_origin

                !(xmin,ymin,zmin) is the l.h.s. corner of the block
                !   -dcellx_grid etc. is added 4/8/2000 T M
                xmin= dcellx_block*(nx-1)
                xmax= xmin+ dcellx_block -dcellx_grid
                ymin= dcelly_block*(ny-1)
                ymax= ymin+ dcelly_block -dcelly_grid
                zmin= dcellz_block*(nz-1)
                zmax= zmin+ dcellz_block -dcellz_grid

                call distsq_blk_atom&
                     (xatom,yatom,zatom,xmin,xmax,ymin,ymax,zmin,zmax,distsq)

                ! if(distsq < rcutsq) then  ! If it is a neighbour block,...
                ! if(distsq < rcutsq+very_small) then  ! If it is a neighbour block,...
                if(distsq < rcutsq(spec)) then  ! If it is a neighbour block,...

                   ! nxmin etc. (those will be used in BtoG trans.)
                   if(naba_blk%no_naba_blk(inp)==0) then
                      naba_blk%nxmin(inp)=nx;naba_blk%nxmax(inp)=nx
                      naba_blk%nymin(inp)=ny;naba_blk%nymax(inp)=ny
                      naba_blk%nzmin(inp)=nz;naba_blk%nzmax(inp)=nz
                   else
                      if(nx<naba_blk%nxmin(inp)) naba_blk%nxmin(inp)=nx
                      if(nx>naba_blk%nxmax(inp)) naba_blk%nxmax(inp)=nx
                      if(ny<naba_blk%nymin(inp)) naba_blk%nymin(inp)=ny
                      if(ny>naba_blk%nymax(inp)) naba_blk%nymax(inp)=ny
                      if(nz<naba_blk%nzmin(inp)) naba_blk%nzmin(inp)=nz
                      if(nz>naba_blk%nzmax(inp)) naba_blk%nzmax(inp)=nz
                   endif

                   ind_block= BCS_blocks%lab_cell(iblock)  !CC in a sim. cell
                   if(ind_block > blocks%mx_gcell) call cq_abort(' ERROR: ind_block in get_naba_BCSblk',ind_block,blocks%mx_gcell)
                   nnd_rem= blocks%i_cc2node(ind_block)  !Domain  responsible node
                   if(nnd_rem > numprocs) call cq_abort(' ERROR : nnd_rem in get_naba_BCSblk',nnd_rem,numprocs)
                   naba_blk%no_naba_blk(inp)=naba_blk%no_naba_blk(inp)+1

                   !checking whether this node is new or not
                    flag_new=.false.
                    if(comm_naba_blocks_of_atoms%no_recv_node(inp) == 0) then
                             flag_new=.true. 
                    elseif(nnd_rem /= &
                        comm_naba_blocks_of_atoms%list_recv_node(comm_naba_blocks_of_atoms%no_recv_node(inp),inp) )  then
                             flag_new=.true. 
                    endif
                   !if new
                   if(flag_new) then
                      comm_naba_blocks_of_atoms%no_recv_node(inp)=comm_naba_blocks_of_atoms%no_recv_node(inp)+1
                      comm_naba_blocks_of_atoms%list_recv_node(comm_naba_blocks_of_atoms%no_recv_node(inp),inp)=nnd_rem
                      !write(io_lun,*) comm_naba_blocks_of_atoms%no_recv_node(inp),'-th node for iprim = ',inp,&
                      !               ' ind_block, nnd_rem = ',ind_block,nnd_rem
                      ncoverx_add =DCS_parts%ncover_rem(1+(nnd_rem-1)*3)
                      ncovery_add =DCS_parts%ncover_rem(2+(nnd_rem-1)*3)
                      ncoverz_add =DCS_parts%ncover_rem(3+(nnd_rem-1)*3)
                      !(difference is between (1-ncoverx_rem) to (ncoverx_rem-1)
                      ncoverz_rem = 2*ncoverz_add+1
                      ncoveryz_rem= ncoverz_rem*(2*ncovery_add+1)

                   endif
                   comm_naba_blocks_of_atoms%no_naba_blk(comm_naba_blocks_of_atoms%no_recv_node(inp),inp)= &
                        comm_naba_blocks_of_atoms%no_naba_blk(comm_naba_blocks_of_atoms%no_recv_node(inp),inp)+1

                   if(naba_blk%no_naba_blk(inp) > naba_blk%mx_naba_blk) then
                      !(NEW) Error check by ierror
                      !    This change is made to print out the information for
                      !   the parameter (mx_naba_blocks_of_atoms).
                      ierror=ierror+1
                   endif

                   nx=anint((xmin-xmin_p+very_small)/dcellx_part)
                   ny=anint((ymin-ymin_p+very_small)/dcelly_part)
                   nz=anint((zmin-zmin_p+very_small)/dcellz_part)
                   nx=nx+ncoverx_add
                   ny=ny+ncovery_add
                   nz=nz+ncoverz_add

                   if(ierror == 0) then
                      !offset from the primary partition to the partition  
                      ! which includes l.h.s. of the covering block.
                      ! It is calculated by using the remote node's ncovers.
                      naba_blk%offset_naba_blk(naba_blk%no_naba_blk(inp),inp)= &
                           (nx-1)*ncoveryz_rem+(ny-1)*ncoverz_rem+nz
                      !CC label in a sim. cell. It will be sent.
                      naba_blk%send_naba_blk  (naba_blk%no_naba_blk(inp),inp)=ind_block 
                      !CC label in a cov. set.  It will be used in BtoG transform.
                      naba_blk%list_naba_blk  (naba_blk%no_naba_blk(inp),inp)=BCS_blocks%lab_cover(iblock)
                   endif ! (ierror=0)
                endif ! (distsq < rcutsq)

             enddo ! iblock in covering sets
             if(iprint_index>3.AND.myid==0) write(io_lun,101) inp,naba_blk%nxmin(inp),naba_blk%nxmax(inp),&
                  naba_blk%nymin(inp),naba_blk%nymax(inp),&
                  naba_blk%nzmin(inp),naba_blk%nzmax(inp)
101          format('inp =',i4,' nxmin,max= ',2i6,' ny ',2i6,' nz ',2i6)

          enddo ! ni (atoms in the primary sets of partitions)

       endif ! Are there atoms ?
    enddo ! np (primary sets of partitions)

    !--------- In case of (ierror >0 ) ----------------
    if(ierror /= 0) then
       inp=0
       do np=1,bundle%groups_on_node  ! primary partitions in bundle
          write(io_lun,*) ' np & atoms in this partition = ',np, bundle%nm_nodgroup(np)
          do ni=1,bundle%nm_nodgroup(np)  ! number of atoms in the partition
             inp=inp+1
             write(io_lun,*) &
                  'inp:naba_blk%no_naba_blk(inp):mx_naba_blk = ', &
                  inp,naba_blk%no_naba_blk(inp),naba_blk%mx_naba_blk
          enddo
       enddo
       call cq_abort('get_naba_BCSblk ierror: ',ierror)
    endif
    return
  end subroutine get_naba_BCSblk
!!***

!!****f* set_blipgrid_module/get_naba_BCSblk_max *
!!
!!  NAME 
!!   get_naba_BCSblk_max
!!  USAGE
!! 
!!  PURPOSE
!!   Finds maxima for BCSblk
!! 
!!  INPUTS
!! 
!!  USES
!! 
!!  AUTHOR
!!   T.Miyazaki/D.R.Bowler
!!  CREATION DATE
!!   2006/09/05
!!  MODIFICATION HISTORY
!!   2016/07/06 17:30 nakata
!!    Renamed max_recv_node_BG -> max_recv_node_BtoG
!!   2016/07/14 16:30 nakata
!!    Renamed max_naba_blk_supp -> max_naba_blocks_of_atoms
!!
!!  SOURCE
!!
  subroutine get_naba_BCSblk_max(rcut,max_naba_blocks_of_atoms, max_recv_node_BtoG)

    use datatypes
    use numbers,        ONLY: very_small, pi,three,four
    use group_module,   ONLY: blocks,parts
    use primary_module, ONLY: bundle
    use cover_module,   ONLY: BCS_blocks,DCS_parts
    use global_module,  ONLY: rcellx,rcelly,rcellz,numprocs
    use block_module,   ONLY: nx_in_block,ny_in_block,nz_in_block
    use GenComms, ONLY: cq_abort
    use species_module, ONLY: n_species

    implicit none

    !Dummy Arguments
    real(double),intent(in),dimension(n_species) :: rcut
    integer :: max_naba_blocks_of_atoms, max_recv_node_BtoG

    !Local variables
    real(double),dimension(n_species) :: rcutsq
    real(double) :: dcellx_part,dcelly_part,dcellz_part
    real(double) :: dcellx_block,dcelly_block,dcellz_block
    real(double) :: dcellx_grid,dcelly_grid,dcellz_grid ! 4/Aug/2000 TM
    real(double) :: xmin_p,ymin_p,zmin_p
    real(double) :: xmin,xmax,ymin,ymax,zmin,zmax
    real(double) :: xatom,yatom,zatom,distsq
    integer :: inp,ncoverz,ncoveryz,np,nx,ny,nz,ni,iblock,ind_block,nnd_rem
    integer :: no_recv_node, no_naba_blk, ii, spec
    integer, dimension(numprocs) :: list_recv_node
    logical :: flag_new

    do ii=1,n_species
       rcutsq(ii)=rcut(ii)*rcut(ii)
    end do
    inp=0 
    max_naba_blocks_of_atoms = 0
    max_recv_node_BtoG = 0
    ncoverz=BCS_blocks%ncoverz
    ncoveryz=BCS_blocks%ncovery*BCS_blocks%ncoverz

    dcellx_part=rcellx/parts%ngcellx ; dcellx_block=rcellx/blocks%ngcellx
    dcelly_part=rcelly/parts%ngcelly ; dcelly_block=rcelly/blocks%ngcelly
    dcellz_part=rcellz/parts%ngcellz ; dcellz_block=rcellz/blocks%ngcellz

    dcellx_grid=dcellx_block/nx_in_block
    dcelly_grid=dcelly_block/ny_in_block
    dcellz_grid=dcellz_block/nz_in_block

    do np=1,bundle%groups_on_node  ! primary partitions in bundle
       if(bundle%nm_nodgroup(np) > 0) then  ! Are there atoms?
          nx=bundle%idisp_primx(np)+bundle%nx_origin
          ny=bundle%idisp_primy(np)+bundle%ny_origin
          nz=bundle%idisp_primz(np)+bundle%nz_origin
          do ni=1,bundle%nm_nodgroup(np)  ! number of atoms in the partition
             inp=inp+1                      ! inp is a primary number of the atom
             if(inp > bundle%mx_iprim) call cq_abort(' Error : iprim in get_naba_BCSblk ',inp, bundle%mx_iprim)
             xatom=bundle%xprim(inp)
             yatom=bundle%yprim(inp)
             zatom=bundle%zprim(inp)
             spec = bundle%species(inp)
             no_naba_blk = 0
             no_recv_node = 0
             list_recv_node = 0
             do iblock=1,BCS_blocks%ng_cover             !NOPG order (cover blks)
                ind_block=BCS_blocks%lab_cover(iblock)-1 !(CC in a cov. set)-1

                nx=ind_block/ncoveryz
                ny=(ind_block-nx*ncoveryz)/ncoverz
                nz=ind_block-nx*ncoveryz-ny*ncoverz
                nx=nx-BCS_blocks%nspanlx+BCS_blocks%nx_origin
                ny=ny-BCS_blocks%nspanly+BCS_blocks%ny_origin
                nz=nz-BCS_blocks%nspanlz+BCS_blocks%nz_origin

                !(xmin,ymin,zmin) is the l.h.s. corner of the block
                !   -dcellx_grid etc. is added 4/8/2000 T M
                xmin= dcellx_block*(nx-1)
                xmax= xmin+ dcellx_block -dcellx_grid
                ymin= dcelly_block*(ny-1)
                ymax= ymin+ dcelly_block -dcelly_grid
                zmin= dcellz_block*(nz-1)
                zmax= zmin+ dcellz_block -dcellz_grid
                call distsq_blk_atom&
                     (xatom,yatom,zatom,xmin,xmax,ymin,ymax,zmin,zmax,distsq)
                if(distsq < rcutsq(spec)) then  ! If it is a neighbour block,...

                   ind_block= BCS_blocks%lab_cell(iblock)  !CC in a sim. cell
                   if(ind_block > blocks%mx_gcell) call cq_abort(' ERROR: ind_block in get_naba_BCSblk',ind_block,blocks%mx_gcell)
                   nnd_rem= blocks%i_cc2node(ind_block)  !Domain  responsible node
                   if(nnd_rem > numprocs) call cq_abort(' ERROR : nnd_rem in get_naba_BCSblk',nnd_rem,numprocs)
                   no_naba_blk=no_naba_blk+1

                   !checking whether this node is new or not
                   flag_new=.false.
                   if(no_recv_node == 0) then
                      flag_new=.true. 
                   else
                      flag_new = .true.
                      do ii = 1,no_recv_node
                         if(nnd_rem==list_recv_node(ii) ) flag_new=.false. 
                      end do
                   endif
                   !if new
                   if(flag_new) then
                      no_recv_node=no_recv_node+1
                      list_recv_node(no_recv_node)=nnd_rem
                   endif
                endif ! (distsq < rcutsq)
             enddo ! iblock in covering sets
             if(no_naba_blk>max_naba_blocks_of_atoms) max_naba_blocks_of_atoms = no_naba_blk
             if(no_recv_node>max_recv_node_BtoG) max_recv_node_BtoG = no_recv_node
          enddo ! ni (atoms in the primary sets of partitions)
       endif ! Are there atoms ?
    enddo ! np (primary sets of partitions)
    
    return
  end subroutine get_naba_BCSblk_max
!!***

!!****f* set_blipgrid_module/distsq_blk_atom *
!!
!!  NAME 
!!   distsq_blk_atom
!!  USAGE
!! 
!!  PURPOSE
!!   sbrt distsq_blk_atom
!!    calculates the square of distance between an atom and a block
!!    We assume a block is orthogonal and it is specified
!!    the minimum and maximum of three direction
!!  (xmin,xmax),(ymin,ymax),(zmin,zmax)
!!    coordinate of atom (xatom,yatom,zatom)
!! 
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!! 
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine distsq_blk_atom &
       (xatom,yatom,zatom,xmin,xmax,ymin,ymax,zmin,zmax,distsq)

    use datatypes
    use numbers, ONLY: ZERO
    implicit none

    real(double),intent(in)::xatom,yatom,zatom
    real(double),intent(in)::xmin,xmax,ymin,ymax,zmin,zmax
    real(double),intent(out)::distsq

    real(double) :: dx,dy,dz

    if(xatom > xmax) then
       dx=xatom-xmax
    elseif(xatom < xmin) then
       dx=xmin-xatom
    else
       dx=ZERO
    endif

    if(yatom > ymax) then
       dy=yatom-ymax
    elseif(yatom < ymin) then
       dy=ymin-yatom
    else
       dy=ZERO
    endif

    if(zatom > zmax) then
       dz=zatom-zmax
    elseif(zatom < zmin) then
       dz=zmin-zatom
    else
       dz=ZERO
    endif

    distsq=dx*dx+dy*dy+dz*dz

    return

  end subroutine distsq_blk_atom
!!***

!!****f* set_blipgrid_module/get_naba_DCSprt *
!!
!!  NAME 
!!   get_naba_DCSprt
!!  USAGE
!! 
!!  PURPOSE
!!   finds neighbour atoms of primary blocks
!!   from DCS_prt covering sets
!!   and makes the list of halo atoms
!! 
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!! 
!!  MODIFICATION HISTORY
!!   2006/07/06 08:27 dave
!!    Changed to use cq_abort
!!  SOURCE
!!
  subroutine get_naba_DCSprt(rcut,naba_set,halo_set)

    use datatypes
    use global_module,  ONLY: rcellx,rcelly,rcellz,numprocs, id_glob, species_glob, sf, nlpf, paof
    use species_module, ONLY: nsf_species, nlpf_species, npao_species, n_species
    use group_module,   ONLY: parts,blocks
    use primary_module, ONLY: domain
    use cover_module,   ONLY: DCS_parts
    use block_module,   ONLY: nx_in_block,ny_in_block,nz_in_block
    use GenComms, ONLY: cq_abort

    implicit none

    real(double),intent(in), dimension(n_species) :: rcut
    type(naba_atm_of_blk),intent(inout) :: naba_set
    type(halo_atm_of_blk),intent(inout) :: halo_set
    real(double) :: xmin,xmax,ymin,ymax,zmin,zmax
    real(double) :: xatom,yatom,zatom,distsq
    real(double), dimension(n_species) :: rcutsq
    real(double) :: dcellx_block,dcelly_block,dcellz_block
    real(double) :: dcellx_grid ,dcelly_grid ,dcellz_grid 
    integer      :: ia,icover,iprim_blk,np,jpart,ni,nnd_old,nnd_rem
    integer      :: ind_part, iorb_alp_i_iblk, iorb_alp_i, norb, spec
    integer      :: irc,ierr

    ! halo_set%naba_atm => naba_set
    do ia = 1,n_species
       rcutsq(ia)=rcut(ia)*rcut(ia)
    end do
    dcellx_block=rcellx/blocks%ngcellx
    dcelly_block=rcelly/blocks%ngcelly
    dcellz_block=rcellz/blocks%ngcellz

    dcellx_grid=dcellx_block/nx_in_block
    dcelly_grid=dcelly_block/ny_in_block
    dcellz_grid=dcellz_block/nz_in_block

    naba_set%no_of_part=0
    naba_set%no_of_atom=0 
    naba_set%ibegin_part(1,1:domain%groups_on_node)=1
    naba_set%no_atom_on_part=0
    
    halo_set%ihalo(:)=0

    !TM VARNSF : START
    naba_set%no_of_orb(:)=0                ! total number of orbitals of naba atoms for this domain
    naba_set%ibegin_blk_orb(:)=0        ! initial position of (alpha, i, iblk) for the block(iblk)
                                        !  in the rank-3 like array
    naba_set%ibeg_orb_atom(:,:)=0     ! initial position of (alpha, i) for the atom(i) in the
                                        !  rank-2 type array (alpha,i) for the block (iblk)
    iorb_alp_i_iblk = 0
    !TM VARNSF : END

    do iprim_blk=1,domain%groups_on_node  ! primary blocks of domain
       !(xmin,ymin,zmin) is the l.h.s. corner of the block
       xmin=(domain%idisp_primx(iprim_blk)+domain%nx_origin-1)*dcellx_block
       ymin=(domain%idisp_primy(iprim_blk)+domain%ny_origin-1)*dcelly_block
       zmin=(domain%idisp_primz(iprim_blk)+domain%nz_origin-1)*dcellz_block
       xmax= xmin+dcellx_block-dcellx_grid
       ymax= ymin+dcelly_block-dcelly_grid
       zmax= zmin+dcellz_block-dcellz_grid

       ia=0                  ! counter of naba atoms for each prim block
       icover=0              ! counter of covering atoms
       jpart=0               ! counter of naba parts for each prim atom

       !TM VARNSF : START
       iorb_alp_i = 0
       naba_set%ibegin_blk_orb(iprim_blk)= iorb_alp_i_iblk + 1
       !TM VARNSF : END

       do np=1,DCS_parts%ng_cover              ! cover set parts (NOPG order)
          jpart=naba_set%no_of_part(iprim_blk)+1     ! index of naba parts
          if(DCS_parts%n_ing_cover(np) > 0) then ! Are there atoms ?
             do ni=1,DCS_parts%n_ing_cover(np)
                icover=icover+1                     ! cover set seq. label of atom
                if(icover /= DCS_parts%icover_ibeg(np)+ni-1 ) &
                     write(io_lun,*) ' ERROR icover np+ni = ', &
                     icover,DCS_parts%icover_ibeg(np)+ni-1

                xatom=DCS_parts%xcover(DCS_parts%icover_ibeg(np)+ni-1)
                yatom=DCS_parts%ycover(DCS_parts%icover_ibeg(np)+ni-1)
                zatom=DCS_parts%zcover(DCS_parts%icover_ibeg(np)+ni-1)

                call distsq_blk_atom &
                     (xatom,yatom,zatom,xmin,xmax,ymin,ymax,zmin,zmax,distsq)
                spec = species_glob( id_glob( parts%icell_beg(DCS_parts%lab_cell(np)) +ni-1 ))
                if(distsq<rcutsq(spec)) then
                !if(distsq < rcutsq) then  ! have found a naba atom
                   ia=ia+1             ! seq. no. of naba atoms for iprim_blk
                   halo_set%ihalo(icover)=1     ! icover-th atom is a halo atom

                   !TM VARNSF : START
                   naba_set%ibeg_orb_atom(ia,iprim_blk)=iorb_alp_i + 1
                   
                   !np in DCS => np_in_sim_cell (through lab_cell)
                   !(np_in_sim_cell, ni) => global_id
                   !norb <= global_id
                   select case(naba_set%function_type)
                   case(sf)
                      norb = nsf_species(spec)
                   case(nlpf)
                      norb = nlpf_species(spec)
                   case(paof)
                      norb = npao_species(spec)
                   case default
                      norb = 0
                   end select
                   iorb_alp_i_iblk = iorb_alp_i_iblk+ norb
                   iorb_alp_i =     iorb_alp_i + norb
                   
                   !TM VARNSF : END
                   if(jpart > naba_set%mx_part) call cq_abort('ERROR: mx_part in get_naba_DCSprt',jpart,naba_set%mx_part)
                   if(ia > naba_set%mx_atom) call cq_abort('ERROR: mx_atom in get_naba_DCSprt',ia,naba_set%mx_atom)
                   naba_set%no_atom_on_part(jpart,iprim_blk)=&
                        naba_set%no_atom_on_part(jpart,iprim_blk)+1
                   naba_set%list_atom(ia,iprim_blk)=ni ! partition seq. label of atom
                   naba_set%list_atom_by_halo(ia,iprim_blk)=icover 
                   !Now in cov. set seq. label -> will be changed to halo seq. label
                endif  ! (distq < rcutsq)
             enddo !  ni (atoms in one of DCS partitions)

             if(jpart > naba_set%mx_part) call cq_abort('ERROR: mx_part in get_naba_DCSprt',jpart,naba_set%mx_part)
             if(naba_set%no_atom_on_part(jpart,iprim_blk)>0) then ! naba part
                naba_set%no_of_part(iprim_blk)      = jpart
                naba_set%list_part(jpart,iprim_blk) = np
                if(jpart < naba_set%mx_part) then
                   naba_set%ibegin_part(jpart+1,iprim_blk)= &
                        naba_set%ibegin_part(jpart,iprim_blk)+ &
                        naba_set%no_atom_on_part(jpart,iprim_blk)
                endif
             endif

          endif ! Are there atoms ?
       enddo ! np (partitions of covering sets)
       naba_set%no_of_atom(iprim_blk)=ia
       naba_set%no_of_orb(iprim_blk)=iorb_alp_i
       if(iprim_blk > 1) then
          naba_set%ibegin_blk(iprim_blk)= &
               naba_set%ibegin_blk(iprim_blk-1)+naba_set%no_of_atom(iprim_blk-1)
       else
          naba_set%ibegin_blk(iprim_blk)=1
       endif
    enddo  ! blocks in domain

!TM VARNSF : START
!    naba_set%no_of_orb=iorb_alp_i_iblk
!TM VARNSF : END

    !
    ! then makes halo atoms, partitions and list of remote nodes
    !
    icover=0
    halo_set%no_of_atom=0
    halo_set%norb=0
    halo_set%no_of_part=0
    halo_set%no_of_node=0
    halo_set%list_of_node(:)=0
    nnd_old=0

    do np=1,DCS_parts%ng_cover          ! Loop over partitions of DCS_parts
       if(DCS_parts%n_ing_cover(np) > 0) then ! Are there atoms in the part?
          ia=0
          do ni=1,DCS_parts%n_ing_cover(np) ! do loop over atoms
             icover=icover+1                 ! index of covering set atoms
             if(icover > halo_set%mx_mcover) call cq_abort(' icover ERROR in get_nabaDCSprt ',icover,halo_set%mx_mcover)
             if(icover /= DCS_parts%icover_ibeg(np)+ni-1) &
                  call cq_abort(' icover ERROR 2 in get_nabaDCSprt ',icover,DCS_parts%icover_ibeg(np)+ni-1)
             if(halo_set%ihalo(icover) == 1) then      ! halo atom
                halo_set%no_of_atom=halo_set%no_of_atom+1
                halo_set%ihalo(icover)=halo_set%no_of_atom
                !TM VARNSF : START
                !SHOULD WE PREPARE # OF ORBITALS for HALO ATOMS?  YES. and we need to know the
                !type of functions. 
                !   It is possible because now I add a member (function_type) in the derived type of
                !   naba_atm_of_blk and halo_set%naba_atm points to naba_atm_of_blk.
                !   itype_func=halo_set%naba_atm%function_type
                spec = species_glob( id_glob( parts%icell_beg(DCS_parts%lab_cell(np)) +ni-1 ))
                norb = 0
                select case(naba_set%function_type)
                case(sf)
                   norb = nsf_species(spec)
                case(nlpf)
                   norb = nlpf_species(spec)
                case(paof)
                   norb = npao_species(spec)
                end select
                if(halo_set%no_of_atom>halo_set%mx_mcover) &
                     call cq_abort("Overflow in get_naba_DCSprt: ",halo_set%no_of_atom,halo_set%mx_mcover)
                halo_set%norb(halo_set%no_of_atom) = norb !<= (np_sim_cell, ni) <= (np, ni)
                !TM VARNSF : END
                ia=ia+1
             endif

          enddo  ! end do loop over atoms in the partition

          if(ia > 0) then ! then halo-partition
             halo_set%no_of_part=halo_set%no_of_part+1
             if(halo_set%no_of_part > halo_set%mx_part) call cq_abort('ERROR in halo_set%mx_part in get_nabaDCSpart',&
                  halo_set%no_of_part,halo_set%mx_part)
             !list of halo parts (NOPG in DCS_parts)
             halo_set%list_of_part(halo_set%no_of_part)=np 
             ind_part=DCS_parts%lab_cell(np)
             nnd_rem=parts%i_cc2node(ind_part)
             if(nnd_rem > numprocs) call cq_abort(' ERROR in nnd_rem from halo_part ', nnd_rem)
             if(halo_set%no_of_part == 1 .OR. nnd_rem /= nnd_old) then !new node
                nnd_old=nnd_rem
                halo_set%no_of_node=halo_set%no_of_node+1
                if(halo_set%no_of_node > halo_set%mx_node) call cq_abort('ERROR! halo_set%mx_node in get_naba_DCSprt ', &
                     halo_set%no_of_node,halo_set%mx_node)
                halo_set%list_of_node(halo_set%no_of_node)=nnd_rem
                !halo_set%ibegin_node(halo_set%no_of_node)=halo_set%no_of_atom-ia+1
             endif  ! if new node

          endif   ! if halo partition

       endif ! Are there atoms in the partition ?
    enddo  ! end do loop over partitions of DCS_parts
    ! Makes (naba_set%list_atom_by_halo) which points naba atms to halo atms
    do iprim_blk=1,domain%groups_on_node         !primary blocks in domain
       if(naba_set%no_of_atom(iprim_blk) > 0) then !Are there naba atoms?
          do ia=1,naba_set%no_of_atom(iprim_blk)     !Loop over naba atoms
             icover=naba_set%list_atom_by_halo(ia,iprim_blk)
             if(icover<1.OR.icover>DCS_parts%mx_mcover ) call cq_abort('Error of icover in get_naba_DCSprt',icover)
             naba_set%list_atom_by_halo(ia,iprim_blk)=halo_set%ihalo(icover)
          enddo
       endif
    enddo
    !
    return
  end subroutine get_naba_DCSprt
!!***

!!****f* set_blipgrid_module/get_naba_DCSprt_max *
!!
!!  NAME 
!!   get_naba_DCSprt_max
!!  USAGE
!! 
!!  PURPOSE
!!   Get maxima for DCSprt
!!  INPUTS
!! 
!!  USES
!! 
!!  AUTHOR
!!   T.Miyazaki/D.R.Bowler
!!  CREATION DATE
!!   July 2000 and 2006/09/06
!!  MODIFICATION HISTORY
!!   2008/05/16 ast
!!    Added timers
!!  SOURCE
!!
  subroutine get_naba_DCSprt_max(rcut,max_naba_prt,max_naba_atm,max_halo_part,max_halo_nodes)

    use datatypes
    use global_module,  ONLY: rcellx,rcelly,rcellz,numprocs, id_glob, species_glob, sf, nlpf, paof
    use species_module, ONLY: nsf_species, nlpf_species, npao_species, n_species
    use group_module,   ONLY: parts,blocks
    use primary_module, ONLY: domain
    use cover_module,   ONLY: DCS_parts
    use block_module,   ONLY: nx_in_block,ny_in_block,nz_in_block
    use GenComms, ONLY: cq_abort

    implicit none

    ! Passed variables
    real(double),intent(in),dimension(n_species) :: rcut
    integer, intent(out) :: max_naba_prt,max_naba_atm,max_halo_part,max_halo_nodes

    ! Local variables
    real(double) :: xmin,xmax,ymin,ymax,zmin,zmax
    real(double) :: xatom,yatom,zatom,distsq
    real(double), dimension(n_species) :: rcutsq
    real(double) :: dcellx_block,dcelly_block,dcellz_block
    real(double) :: dcellx_grid ,dcelly_grid ,dcellz_grid 
    integer      :: ia,icover,iprim_blk,np,jpart,ni,nnd_old,nnd_rem
    integer      :: ind_part, iorb_alp_i_iblk, iorb_alp_i, norb, spec
    integer      :: irc,ierr, no_of_part
    integer, allocatable, dimension(:) :: ihalo
    logical :: atoms

    do ia=1,n_species
       rcutsq(ia)=rcut(ia)*rcut(ia)
    end do
    call start_timer(tmr_std_allocation)
    allocate(ihalo(DCS_parts%mx_mcover),STAT=ierr)
    if(ierr/=0) call cq_abort("Error allocating icover in getDCSprtmax: ",DCS_parts%mx_mcover,ierr)
    call stop_timer(tmr_std_allocation)
    dcellx_block=rcellx/blocks%ngcellx
    dcelly_block=rcelly/blocks%ngcelly
    dcellz_block=rcellz/blocks%ngcellz

    dcellx_grid=dcellx_block/nx_in_block
    dcelly_grid=dcelly_block/ny_in_block
    dcellz_grid=dcellz_block/nz_in_block

    iorb_alp_i_iblk = 0
    !TM VARNSF : END
    max_naba_prt = 0
    max_naba_atm = 0
    max_halo_part = 0
    ihalo = 0
    do iprim_blk=1,domain%groups_on_node  ! primary blocks of domain
       !(xmin,ymin,zmin) is the l.h.s. corner of the block
       xmin=(domain%idisp_primx(iprim_blk)+domain%nx_origin-1)*dcellx_block
       ymin=(domain%idisp_primy(iprim_blk)+domain%ny_origin-1)*dcelly_block
       zmin=(domain%idisp_primz(iprim_blk)+domain%nz_origin-1)*dcellz_block
       xmax= xmin+dcellx_block-dcellx_grid
       ymax= ymin+dcelly_block-dcelly_grid
       zmax= zmin+dcellz_block-dcellz_grid

       ia=0                  ! counter of naba atoms for each prim block
       icover=0              ! counter of covering atoms
       jpart=0               ! counter of naba parts for each prim atom
       no_of_part=0 ! Serves purpose of naba_set%no_of_part(iprim_blk)
       do np=1,DCS_parts%ng_cover              ! cover set parts (NOPG order)
          jpart=no_of_part+1     ! index of naba parts
          atoms = .false.
          if(DCS_parts%n_ing_cover(np) > 0) then ! Are there atoms ?
             do ni=1,DCS_parts%n_ing_cover(np)
                icover=icover+1                     ! cover set seq. label of atom
                if(icover /= DCS_parts%icover_ibeg(np)+ni-1 ) &
                     write(io_lun,*) ' ERROR icover np+ni = ', &
                     icover,DCS_parts%icover_ibeg(np)+ni-1

                xatom=DCS_parts%xcover(DCS_parts%icover_ibeg(np)+ni-1)
                yatom=DCS_parts%ycover(DCS_parts%icover_ibeg(np)+ni-1)
                zatom=DCS_parts%zcover(DCS_parts%icover_ibeg(np)+ni-1)

                call distsq_blk_atom(xatom,yatom,zatom,xmin,xmax,ymin,ymax,zmin,zmax,distsq)
                spec = species_glob( id_glob( parts%icell_beg(DCS_parts%lab_cell(np)) +ni-1 ))

                if(distsq < rcutsq(spec)) then  ! have found a naba atom
                   ia=ia+1          
                   atoms = .true.
                   ihalo(icover) = 1
                endif  ! (distq < rcutsq)
             enddo !  ni (atoms in one of DCS partitions)
             if(atoms) then 
                no_of_part = jpart
             end if
          endif ! Are there atoms ?
       enddo ! np (partitions of covering sets)
       if(ia>max_naba_atm) max_naba_atm = ia
       if(jpart>max_naba_prt) max_naba_prt = jpart
    enddo  ! blocks in domain
    ! Find partitions in halo
    icover=0
    nnd_old = 0
    max_halo_nodes = 0
    do np=1,DCS_parts%ng_cover          ! Loop over partitions of DCS_parts
       if(DCS_parts%n_ing_cover(np) > 0) then ! Are there atoms in the part?
          ia=0
          do ni=1,DCS_parts%n_ing_cover(np) ! do loop over atoms
             icover=icover+1                 ! index of covering set atoms
             if(ihalo(icover) == 1) ia=ia+1
          enddo  ! end do loop over atoms in the partition
          if(ia > 0) then
             max_halo_part = max_halo_part + 1
             ind_part=DCS_parts%lab_cell(np)
             nnd_rem=parts%i_cc2node(ind_part)
             if(nnd_rem > numprocs) call cq_abort(' ERROR in nnd_rem from halo_part ', nnd_rem)
             if(nnd_rem /= nnd_old) then !new node
                nnd_old=nnd_rem
                max_halo_nodes=max_halo_nodes+1
             endif  ! if new node
          end if
       endif ! Are there atoms in the partition ?
    enddo  ! end do loop over partitions of DCS_parts
    call start_timer(tmr_std_allocation)
    deallocate(ihalo,STAT=ierr)
    if(ierr/=0) call cq_abort("Error deallocating icover in getDCSprtmax: ",DCS_parts%mx_mcover,ierr)
    call stop_timer(tmr_std_allocation)
    return
  end subroutine get_naba_DCSprt_max
!!***

!!****f* set_blipgrid_module/make_sendinfo_BtoG *
!!
!!  NAME 
!!   make_sendinfo_BtoG
!!  USAGE
!! 
!!  PURPOSE
!!   prepares information of receiving in Blip-Grid transforms
!! = makes comm_naba_blocks_of_atoms from naba_atoms_of_blocks(supp)
!! 
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   T.Miyazaki
!!  CREATION DATE
!! 
!!  MODIFICATION HISTORY
!!   08:18, 2003/06/11 dave
!!    Added TM's debug statements
!!   11:55, 12/11/2004 dave 
!!    Changed distribute_atom_module to atoms
!!   2016/07/06 17:30 nakata
!!    Renamed subroutine make_sendinfo_BG -> make_sendinfo_BtoG 
!!    Renamed comm_in_BG -> comm_in_BtoG and comBG -> comm_naba_blocks_of_atoms
!!   2016/07/14 16:30 nakata
!!    Renamed naba_supp -> naba_atm_set
!!  SOURCE
!!
  subroutine make_sendinfo_BtoG(myid,naba_atm_set,comm_naba_blocks_of_atoms)

    use datatypes
    use global_module,          ONLY:id_glob
    use atoms, ONLY:atom_number_on_node
    use group_module,           ONLY:parts
    use primary_module,         ONLY:domain
    use cover_module,           ONLY:DCS_parts
    use naba_blk_module,        ONLY:naba_atm_of_blk,comm_in_BtoG
    !08/04/2003 T. Miyazaki
    use global_module,          ONLY: iprint_index
!    use maxima_module,          ONLY: mx_sent_pair_BtoG, mx_send_node_BtoG &
!                                    , mx_recv_node_BtoG
    use GenComms,               ONLY: cq_abort, my_barrier
    implicit none
    type(naba_atm_of_blk),intent(in) ::naba_atm_set
    type(comm_in_BtoG),intent(inout)::comm_naba_blocks_of_atoms
    integer,intent(in) :: myid

    integer :: iprim_blk,ipart,ipair,iprim,inode
    integer :: jpart,ind_part,nnd_rem,ibegin,iend
    integer :: ni,ia,iprim_rem
    integer :: ifind,ii,ind_node,isend,npair
    integer :: irc,ierr
    integer :: mx_sent_pair,mx_send_node,mx_recv_node
    !08/04/2003 T. Miyazaki
    integer, parameter :: iprint_debug = 10
    integer, parameter :: iprint_debug2 = 20
 
    comm_naba_blocks_of_atoms%no_send_node(:)=0
    comm_naba_blocks_of_atoms%no_sent_pairs(:,:)=0
    comm_naba_blocks_of_atoms%list_send_node(:,:)=0

    do iprim_blk=1,domain%groups_on_node          ! Do loop over prim blks
       if(naba_atm_set%no_of_part(iprim_blk) > 0) then ! are there naba atoms?
          do ipart=1,naba_atm_set%no_of_part(iprim_blk)   ! naba parts of prim blk
             jpart=naba_atm_set%list_part(ipart,iprim_blk)
             ind_part=DCS_parts%lab_cell(jpart)
             nnd_rem=parts%i_cc2node(ind_part)

             if(naba_atm_set%no_atom_on_part(ipart,iprim_blk) < 1) then  
                !!   check no_naba_atom
                write(io_lun,*) 'no of atoms in the neighbour partition is 0 ???'
                write(io_lun,*) ' ERROR in make_sendinfo_BtoG &
                     &for iprim_blk,ipart,jpart = ' &
                     ,iprim_blk,ipart,jpart
                call cq_abort('ERROR in naba_atm_set%no_atom_on_part in make_sendinfo_BtoG')
             endif
             ibegin=naba_atm_set%ibegin_part(ipart,iprim_blk)
             iend  =ibegin + naba_atm_set%no_atom_on_part(ipart,iprim_blk)-1

             do ipair=ibegin,iend                        ! naba atms in the part
                ni=naba_atm_set%list_atom(ipair,iprim_blk)
                ia=id_glob(parts%icell_beg(ind_part)+ni-1) !global id
                iprim_rem= atom_number_on_node(ia)

                !(check existed sending nodes for iprim_rem)
                ifind=0
                if(comm_naba_blocks_of_atoms%no_send_node(iprim_rem) > 0) then  
                   do ii=1,comm_naba_blocks_of_atoms%no_send_node(iprim_rem)
                      ind_node=comm_naba_blocks_of_atoms%list_send_node(ii,iprim_rem)
                      if(nnd_rem == ind_node) then
                         ifind=ii
                         exit
                      endif
                   enddo
                endif ! (check  existed sending nodes for iprim_rem)

                if(ifind == 0) then   ! (if the current sending node is new )
                   comm_naba_blocks_of_atoms%no_send_node(iprim_rem)=comm_naba_blocks_of_atoms%no_send_node(iprim_rem)+1
                   if(comm_naba_blocks_of_atoms%no_send_node(iprim_rem) > comm_naba_blocks_of_atoms%mx_send_node) then
                      write(io_lun,*) ' ERROR mx_send_node in make_sendinfo_BtoG', &
                           comm_naba_blocks_of_atoms%no_send_node(iprim_rem), comm_naba_blocks_of_atoms%mx_send_node
                      call cq_abort('ERROR in mx_send_node in make_sendinfo_BtoG')
                   endif
                   ! write(io_lun,*) 'no_send_node = ',comm_naba_blocks_of_atoms%no_send_node(iprim_rem)
                   ! write(io_lun,*) 'B list_send_node', &
                   !  comm_naba_blocks_of_atoms%list_send_node(comm_naba_blocks_of_atoms%no_send_node(iprim_rem),iprim_rem), &
                   !         'no_sent_pairs for this node at present ', &
                   !  comm_naba_blocks_of_atoms%no_sent_pairs(comm_naba_blocks_of_atoms%no_send_node(iprim_rem),iprim_rem) 

                   comm_naba_blocks_of_atoms%list_send_node &
                        (comm_naba_blocks_of_atoms%no_send_node(iprim_rem),iprim_rem)=nnd_rem
                   comm_naba_blocks_of_atoms%no_sent_pairs(comm_naba_blocks_of_atoms%no_send_node(iprim_rem),iprim_rem)= &
                        comm_naba_blocks_of_atoms%no_sent_pairs(comm_naba_blocks_of_atoms%no_send_node(iprim_rem),iprim_rem) + 1

                   !write(io_lun,*) 'A list_send_node', &
                   !    comm_naba_blocks_of_atoms%list_send_node(comm_naba_blocks_of_atoms%no_send_node(iprim_rem),iprim_rem), &
                   !           'no_sent_pairs for this node at present ', &
                   !    comm_naba_blocks_of_atoms%no_sent_pairs(comm_naba_blocks_of_atoms%no_send_node(iprim_rem),iprim_rem) 

                else                  ! (if the current sending node is existed )
                   comm_naba_blocks_of_atoms%no_sent_pairs(ifind,iprim_rem)= &
                        comm_naba_blocks_of_atoms%no_sent_pairs(ifind,iprim_rem) + 1
                endif                 !(the present node is new or old)

                ! write(io_lun,101) myid,&
                !  nnd_rem,ni,ia,iprim_rem, comm_naba_blocks_of_atoms%no_send_node(iprim_rem), &
                !  (comm_naba_blocks_of_atoms%no_sent_pairs(ii,iprim_rem),ii=1,comm_naba_blocks_of_atoms%no_send_node(iprim_rem))
                !101 format(i3,'RemNode',i3,' ni, ia = ',i3,i5,'iprim_rem= ',i4, &
                !             ' no_send_node&pairs (iprim_rem) = ',i3,3x,4i5)


             enddo  ! neighbour atoms in the partition
          enddo   ! neighbour partition of primary blocks
       endif   ! Are there neighbour atoms (partitions) ?
    enddo    ! primary blocks
    !check comm_naba_blocks_of_atoms%mx_pair and comm_naba_blocks_of_atoms%mx_recv_call
!%%!    isend=0;npair=0
!%%!    do iprim=1,comm_naba_blocks_of_atoms%mx_iprim
!%%!       if(comm_naba_blocks_of_atoms%no_send_node(iprim) > 0) then
!%%!          do inode=1,comm_naba_blocks_of_atoms%no_send_node(iprim)
!%%!             isend=isend+1
!%%!             if(isend == 1) then
!%%!                mx_sent_pair=comm_naba_blocks_of_atoms%no_sent_pairs(inode,iprim)
!%%!             elseif(comm_naba_blocks_of_atoms%no_sent_pairs(inode,iprim) > mx_sent_pair) then
!%%!                mx_sent_pair=comm_naba_blocks_of_atoms%no_sent_pairs(inode,iprim)
!%%!             endif
!%%!
!%%!             npair=npair+comm_naba_blocks_of_atoms%no_sent_pairs(inode,iprim)
!%%!             if(comm_naba_blocks_of_atoms%no_sent_pairs(inode,iprim) > comm_naba_blocks_of_atoms%mx_pair) then
!%%!                write(io_lun,*) ' ERROR mx_pair in make_sendinfo_BtoG'
!%%!                call cq_abort('ERROR mx_pair in make_sendinfo_BtoG')
!%%!             endif
!%%!          enddo
!%%!       endif !(comm_naba_blocks_of_atoms%no_send_node(iprim) > 0) 
!%%!    enddo
!%%!    do iprim=1,comm_naba_blocks_of_atoms%mx_iprim
!%%!       if(iprim==1) then
!%%!          mx_send_node=comm_naba_blocks_of_atoms%no_send_node(iprim)
!%%!       elseif(comm_naba_blocks_of_atoms%no_send_node(iprim) > mx_send_node) then
!%%!          mx_send_node=comm_naba_blocks_of_atoms%no_send_node(iprim)
!%%!       endif
!%%!    enddo
!%%!    do iprim=1,comm_naba_blocks_of_atoms%mx_iprim
!%%!       if(iprim==1) then
!%%!          mx_recv_node=comm_naba_blocks_of_atoms%no_recv_node(iprim)
!%%!       elseif(comm_naba_blocks_of_atoms%no_recv_node(iprim) > mx_recv_node) then
!%%!          mx_recv_node=comm_naba_blocks_of_atoms%no_recv_node(iprim)
!%%!       endif
!%%!    enddo
!%%!
!%%!    if(iprint > iprint_debug ) then
!%%!       write(io_lun,*) ' @@ mx_sent_pair_BtoG should be ',mx_sent_pair,&
!%%!            ' for node ',myid+1
!%%!       write(io_lun,*) ' @@ mx_send_node_BtoG should be ',mx_send_node,&
!%%!            ' for node ',myid+1
!%%!       write(io_lun,*) ' @@ mx_recv_node_BtoG should be ',mx_recv_node,&
!%%!            ' for node ',myid+1
!%%!    endif
!%%!    ! To check mx_sent_pair_BtoG etc.
!%%!    ierr = 0
!%%!    if(mx_sent_pair > mx_sent_pair_BtoG) then
!%%!       ierr = ierr + 1   
!%%!       write(io_lun,*) 'ERROR! mx_sent_pair_BtoG must be larger than ',&
!%%!            mx_sent_pair,' present value = ', mx_sent_pair_BtoG
!%%!    endif
!%%!    if(mx_send_node > mx_send_node_BtoG) then
!%%!       ierr = ierr + 2
!%%!       write(io_lun,*) 'ERROR! mx_send_node_BtoG must be larger than ',&
!%%!            mx_send_node,' present value = ', mx_send_node_BtoG
!%%!    endif
!%%!    if(mx_recv_node > mx_recv_node_BtoG) then
!%%!       ierr = ierr + 4
!%%!       write(io_lun,*) 'ERROR! mx_recv_node_BtoG must be larger than ',&
!%%!            mx_recv_node,' present value = ', mx_recv_node_BtoG
!%%!    endif
!%%!    call my_barrier()
!%%!    if(ierr /= 0) then
!%%!       call cq_abort('Error in setting mx_.._BtoG problem',ierr)
!%%!    endif
!%%!
!%%!    if(iprint > iprint_debug2) then
!%%!       write(io_lun,*) ' Node ', myid+1,'  isend,npair = ',isend,npair
!%%!       write(io_lun,*) 'COMBtoGMXIPRIM= ',comm_naba_blocks_of_atoms%mx_iprim
!%%!       do iprim=1,comm_naba_blocks_of_atoms%mx_iprim
!%%!          write(io_lun,*) myid+1,' IPRIM = ',iprim,' No of remote nodes = ', &
!%%!               comm_naba_blocks_of_atoms%no_send_node(iprim),' LIST SEND= ', &
!%%!               (comm_naba_blocks_of_atoms%list_send_node(inode,iprim),inode=1,comm_naba_blocks_of_atoms%no_send_node(iprim))&
!%%!               ,' LIST RECV= ',&
!%%!               (comm_naba_blocks_of_atoms%list_recv_node(inode,iprim),inode=1,comm_naba_blocks_of_atoms%no_recv_node(iprim))
!%%!       enddo
!%%!    endif

    return 
  end subroutine make_sendinfo_BtoG
!!***

!!****f* set_blipgrid_module/make_sendinfo_BtoG_max *
!!
!!  NAME 
!!   make_sendinfo_BtoG_max
!!  USAGE
!! 
!!  PURPOSE
!!   Finds maxima for comm_naba_blocks_of_atoms
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   T.Miyazaki/D.R.Bowler
!!  CREATION DATE
!!   June 2000 and 2006/09/06
!!  MODIFICATION HISTORY
!!   2008/05/16 ast
!!    Added timers
!!   2016/07/06 17:30 nakata
!!    Renamed subroutine make_sendinfo_BG_max -> make_sendinfo_BtoG_max
!!    Renamed comm_in_BG -> comm_in_BtoG
!!   2016/07/14 16:30 nakata
!!    Renamed naba_supp -> naba_atm_set
!!  SOURCE
!!
  subroutine make_sendinfo_BtoG_max(myid,naba_atm_set,max_send_node,max_sent_pairs,max_recv_call)

    use datatypes
    use global_module,          ONLY:id_glob
    use atoms, ONLY:atom_number_on_node
    use group_module,           ONLY:parts
    use primary_module,         ONLY:domain,bundle
    use cover_module,           ONLY:DCS_parts
    use naba_blk_module,        ONLY:naba_atm_of_blk,comm_in_BtoG
    !08/04/2003 T. Miyazaki
    use global_module,          ONLY: iprint_index, numprocs
!    use maxima_module,          ONLY: mx_sent_pair_BtoG, mx_send_node_BtoG &
!                                    , mx_recv_node_BtoG
    use GenComms,               ONLY: cq_abort, my_barrier
    implicit none
    type(naba_atm_of_blk),intent(in) ::naba_atm_set
    integer, intent(out) :: max_send_node,max_sent_pairs,max_recv_call
    integer,intent(in) :: myid

    integer :: iprim_blk,ipart,ipair,iprim,inode
    integer :: jpart,ind_part,nnd_rem,ibegin,iend
    integer :: ni,ia,iprim_rem
    integer :: ifind,ii,ind_node,isend,npair
    integer :: irc,ierr
    integer :: mx_sent_pair,mx_send_node,mx_recv_node
    integer, allocatable, dimension(:) :: no_send_node
    integer, allocatable, dimension(:,:) :: list_send_node, no_sent_pairs
 
    max_recv_call = 0
    max_send_node = 0
    max_sent_pairs = 0
    call start_timer(tmr_std_allocation)
    allocate(no_send_node(bundle%mx_iprim),STAT=ierr)
    if(ierr/=0) call cq_abort("Error allocating no_send_node in BtoGmax: ",bundle%mx_iprim,ierr)
    allocate(list_send_node(numprocs,bundle%mx_iprim),STAT=ierr)
    if(ierr/=0) call cq_abort("Error allocating list_send_node in BtoGmax: ",bundle%mx_iprim,numprocs)
    allocate(no_sent_pairs(numprocs,bundle%mx_iprim),STAT=ierr)
    if(ierr/=0) call cq_abort("Error allocating no_sent_pairs in BtoGmax: ",bundle%mx_iprim,numprocs)
    call stop_timer(tmr_std_allocation)
    no_send_node = 0
    list_send_node = 0
    no_sent_pairs=0
    do iprim_blk=1,domain%groups_on_node          ! Do loop over prim blks
       if(naba_atm_set%no_of_part(iprim_blk) > 0) then ! are there naba atoms?
          do ipart=1,naba_atm_set%no_of_part(iprim_blk)   ! naba parts of prim blk
             jpart=naba_atm_set%list_part(ipart,iprim_blk)
             ind_part=DCS_parts%lab_cell(jpart)
             nnd_rem=parts%i_cc2node(ind_part)
             if(naba_atm_set%no_atom_on_part(ipart,iprim_blk) < 1) then  
                !!   check no_naba_atom
                write(io_lun,*) 'no of atoms in the neighbour partition is 0 ???'
                write(io_lun,*) ' ERROR in make_sendinfo_BtoG for iprim_blk,ipart,jpart = ' ,iprim_blk,ipart,jpart
                call cq_abort('ERROR in naba_atm_set%no_atom_on_part in make_sendinfo_BtoG')
             endif
             ibegin=naba_atm_set%ibegin_part(ipart,iprim_blk)
             iend  =ibegin + naba_atm_set%no_atom_on_part(ipart,iprim_blk)-1
             do ipair=ibegin,iend                        ! naba atms in the part
                ni=naba_atm_set%list_atom(ipair,iprim_blk)
                ia=id_glob(parts%icell_beg(ind_part)+ni-1) !global id
                iprim_rem= atom_number_on_node(ia)
                ifind=0
                if(no_send_node(iprim_rem) > 0) then  
                   do ii=1,no_send_node(iprim_rem)
                      ind_node=list_send_node(ii,iprim_rem)
                      if(nnd_rem == ind_node) then
                         ifind=ii
                         exit
                      endif
                   enddo
                endif ! (check  existed sending nodes for iprim_rem)
                if(ifind == 0) then   ! (if the current sending node is new )
                   no_send_node(iprim_rem)=no_send_node(iprim_rem)+1
                   list_send_node(no_send_node(iprim_rem),iprim_rem)=nnd_rem
                   no_sent_pairs(no_send_node(iprim_rem),iprim_rem)= &
                        no_sent_pairs(no_send_node(iprim_rem),iprim_rem) + 1
                else                  ! (if the current sending node is existed )
                   no_sent_pairs(ifind,iprim_rem)= no_sent_pairs(ifind,iprim_rem) + 1
                endif                 !(the present node is new or old)
             enddo  ! neighbour atoms in the partition
          enddo   ! neighbour partition of primary blocks
       endif   ! Are there neighbour atoms (partitions) ?
    enddo    ! primary blocks
    do iprim=1,bundle%mx_iprim
       if(no_send_node(iprim)>max_send_node) max_send_node = no_send_node(iprim)
       if(no_send_node(iprim) > 0) then
          do inode=1,no_send_node(iprim)
             if(no_sent_pairs(inode,iprim)>max_sent_pairs) max_sent_pairs = no_sent_pairs(inode,iprim)
          end do
       endif ! (no_send_node(iprim) > 0) then
       max_recv_call = max_recv_call + no_send_node(iprim)
    end do
    call start_timer(tmr_std_allocation)
    deallocate(no_sent_pairs, list_send_node, no_send_node,STAT=ierr)
    if(ierr/=0) call cq_abort("Error deallocating no_sent_pairs in BtoGmax: ",bundle%mx_iprim,numprocs)
    call stop_timer(tmr_std_allocation)
    return 
  end subroutine make_sendinfo_BtoG_max
!!***

!!****f* set_blipgrid_module/make_table_BtoG *
!!
!!  NAME 
!!   make_table_BtoG
!!  USAGE
!! 
!!  PURPOSE
!!   makes table showing where to put sent BtoG transformed functions
!!   Its main task (done in subroutine make_table) is to calculate 
!!
!!      comm_naba_blocks_of_atoms%table_blk (ipair,isend) : seq. no. of primary block
!!      comm_naba_blocks_of_atoms%table_pair(ipair,isend) : seq. no. of naba atom for 
!!                                      the above primary block
!!
!!    Here, ipair is the seq. no. within each receive call and
!!    isend is the seq. no. of my receiving calls.
!! 
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   T. Miyazaki
!!  CREATION DATE
!!   2000/09/08
!!  MODIFICATION HISTORY
!!   2006/07/06 08:30 dave
!!    Tidying
!!   2016/07/06 17:30 nakata
!!    Renamed subroutine make_table_BG -> make_table_BtoG
!!    Renamed comm_in_BG -> comm_in_BtoG and comBG to comm_naba_blocks_of_atoms
!!  SOURCE
!!
  subroutine make_table_BtoG(myid,comm_naba_blocks_of_atoms)

    use mpi, ONLY: MPI_STATUS_SIZE
    use datatypes
    use maxima_module, ONLY: maxnsf
    use primary_module,  ONLY:bundle
    use block_module, ONLY: n_pts_in_block
    use GenComms, ONLY: my_barrier, cq_abort

    implicit none

    !Dummy Arguments
    integer,intent(in) :: myid
    !type(naba_blk_of_atm),intent(in) :: naba_blk
    !type(naba_atm_of_blk),intent(in) :: naba_atoms_of_blocks(supp)
    type(comm_in_BtoG),intent(inout)   :: comm_naba_blocks_of_atoms
    !Local variables           
    integer,parameter :: tag=1  
    !integer,save :: iprim
    integer :: iprim, j
    integer :: ierr,irc,istart_myid
    integer :: mynode
    integer :: isend, stat
    integer :: mx_send,mx_recv
    integer :: nsend_req(comm_naba_blocks_of_atoms%mx_recv_node)
    integer, dimension(MPI_STATUS_SIZE) :: mpi_stat
    
    integer :: ii
    mynode=myid+1
    isend=0           ! index of receiving data from sending nodes

    call my_barrier()
    do iprim=1,bundle%mx_iprim
       ! set comm_naba_blocks_of_atoms%ibegin_recv_call(iprim)
       if(iprim == 1) then
          comm_naba_blocks_of_atoms%ibeg_recv_call(iprim)=1
       else
          comm_naba_blocks_of_atoms%ibeg_recv_call(iprim)=isend+1
       endif
       if(iprim <=  bundle%n_prim) then
          call send_array_BtoG(iprim,mynode,istart_myid,nsend_req)
       else
          istart_myid=0
       endif
       call recv_array_BtoG(iprim,mynode,istart_myid,isend)
       do j = 1,comm_naba_blocks_of_atoms%no_recv_node(iprim)
          if(comm_naba_blocks_of_atoms%list_recv_node(j,iprim)/=mynode) call MPI_Wait(nsend_req(j),mpi_stat,ierr)
       end do
       call my_barrier()
       ! we need this BARRIER as the send_array of a node which finish
       ! receiving MPI calls might be changed before other nodes have
       ! not finished receiving this array.
       if(allocated(isend_array)) then
          call start_timer(tmr_std_allocation)
          deallocate(isend_array,STAT=stat)
          if(stat/=0) call cq_abort("Error deallocating isend_array in send_array_BtoG: ",stat)
          call stop_timer(tmr_std_allocation)
       end if
    enddo
    return
  end subroutine make_table_BtoG
!!***

  !-------------------------------------------------------------------
  !sbrt send_array_BtoG
  ! send no_naba_blks for each node,
  ! send_naba_blk and offset_naba_blk
  !-------------------------------------------------------------------
  subroutine send_array_BtoG(iprim,mynode,istart_myid,nsend_req)

    !Modules and Dummy Arguments --
    use mpi, ONLY: MPI_INTEGER, MPI_COMM_WORLD
    use GenComms, ONLY: cq_abort

    implicit none

    integer,intent(in)  :: iprim,mynode
    integer,intent(out) :: istart_myid
    integer,intent(out) :: nsend_req(comm_naba_blocks_of_atoms%mx_recv_node)
    !Local variables -- 
    integer :: nnodes,istart,inode,nnd,ibegin
    integer :: iend,nsize,send_size, tmpsize, stat
    integer :: ierr(comm_naba_blocks_of_atoms%mx_recv_node)
    integer,parameter :: tag=1

    nnodes=comm_naba_blocks_of_atoms%no_recv_node(iprim)

    istart=1;istart_myid=0
    if(nnodes > 0) then
       tmpsize = 0
       do inode=1,nnodes
          nsize =comm_naba_blocks_of_atoms%no_naba_blk(inode,iprim)
          tmpsize = tmpsize+2*nsize+1
       enddo
       if(allocated(isend_array)) then 
          write(io_lun,*) mynode,' send_array alloc: ',size(isend_array)
          call cq_abort('Problem with send_array !')
       end if
       call start_timer(tmr_std_allocation)
       allocate(isend_array(tmpsize),STAT=stat)
       if(stat/=0) call cq_abort("Error allocating isend_array in send_array_BtoG: ",tmpsize)
       call stop_timer(tmr_std_allocation)
       do inode=1,nnodes
          nnd=comm_naba_blocks_of_atoms%list_recv_node(inode,iprim)
          if(inode == 1) then
             ibegin=1
          else      
             ibegin=iend+1
          endif
          nsize =comm_naba_blocks_of_atoms%no_naba_blk(inode,iprim)
          iend  =ibegin+nsize-1

          isend_array(istart)=nsize
          isend_array(istart+        1:istart+  nsize)=&
               naba_blocks_of_atoms%send_naba_blk  (ibegin:iend,iprim)
          isend_array(istart+  nsize+1:istart+2*nsize)=&
               naba_blocks_of_atoms%offset_naba_blk(ibegin:iend,iprim)
          send_size=2*nsize+1

          if(nnd /= mynode) then
             call MPI_issend(isend_array(istart),send_size,MPI_INTEGER,nnd-1, &
                  tag,MPI_COMM_WORLD,nsend_req(inode),ierr(inode))
          else ! nnd==mynode
             istart_myid=istart
          endif
          istart=istart+2*nsize+1
       enddo
    endif
    return
  end subroutine send_array_BtoG

  !-------------------------------------------------------------------
  !sbrt recv_array_BtoG
  !  Domain responsible node receives the data from sending nodes
  ! and makes the table for storing the received blip-grid 
  ! transformed support functions. (See subroutine make_table)
  ! 
  !  We have to call (make_table) from the do-loop
  ! in this subroutine to reuse the array irecv_array.
  !-------------------------------------------------------------------
  subroutine recv_array_BtoG(iprim,mynode,istart_myid,isend)

    use global_module, ONLY: numprocs
    use primary_module, ONLY: bundle
    use mpi, ONLY: MPI_STATUS_SIZE, MPI_INTEGER, MPI_COMM_WORLD
    use GenComms, ONLY: cq_abort

    implicit none

    !Dummy Arguments
    integer,intent(in)   ::iprim,mynode,istart_myid
    integer,intent(inout)::isend
    !Local variables
    integer::nnodes,nnd_rem,inode,npair,recv_size,stat
    integer,parameter :: tag=1
    integer :: nrecv_stat(MPI_STATUS_SIZE,comm_naba_blocks_of_atoms%mx_send_node)
    integer :: ierr(comm_naba_blocks_of_atoms%mx_send_node)
    integer:: irc

    nnodes=comm_naba_blocks_of_atoms%no_send_node(iprim)
    if(nnodes > numprocs .or. nnodes < 0) call cq_abort(' no of nodes ERROR ', &
                                                        comm_naba_blocks_of_atoms%no_send_node(iprim), numprocs )
    !write(mynode+19,*) 'Node ',mynode,' Remote Nodes ',nnodes,' iprim = ',iprim
    !do inode=1,nnodes
    ! write(mynode+19,*) '   inode, nnd_rem = ',inode,comm_naba_blocks_of_atoms%list_send_node(inode,iprim)
    !enddo

    if(nnodes > 0) then
       do inode=1,nnodes

          !STRANGE AGAIN !! ----  FIXED
          if(inode > nnodes) call cq_abort('WARNING!!   in recv_array_BtoG', inode,nnodes)
          nnd_rem=comm_naba_blocks_of_atoms%list_send_node(inode,iprim)
          recv_size=comm_naba_blocks_of_atoms%no_sent_pairs(inode,iprim)*2+1
          call start_timer(tmr_std_allocation)
          allocate(irecv_array(recv_size),STAT=stat)
          if(stat/=0) call cq_abort("Error allocating irecv_array in recv_array_BtoG: ",recv_size)
          call stop_timer(tmr_std_allocation)
          if(nnd_rem /= mynode) then
             call MPI_recv(irecv_array,recv_size,MPI_INTEGER,nnd_rem-1, &
                  tag,MPI_COMM_WORLD,nrecv_stat(:,inode),ierr(inode))
             if(ierr(inode) /= 0)  call cq_abort('MPI_recv: ierr = ',ierr(inode))
             npair=irecv_array(1)
             naba_blk_rem => irecv_array(        2:  npair+1)
             offset_rem => irecv_array(  npair+2:2*npair+1)

          else

             !write(io_lun,*) ' Recv Node ',mynode,' Istart: ',istart_myid
             npair=isend_array(istart_myid)
             naba_blk_rem => isend_array(istart_myid+        1:istart_myid+npair)
             offset_rem => isend_array(istart_myid+  npair+1:istart_myid+2*npair)

          endif

          !Check npair between sending and receiving nodes
          if(npair /= comm_naba_blocks_of_atoms%no_sent_pairs(inode,iprim)) call cq_abort('pair problem in recv_array_BtoG: ', &
               npair,comm_naba_blocks_of_atoms%no_sent_pairs(inode,iprim))
          isend=isend+1
          if(isend > comm_naba_blocks_of_atoms%mx_recv_call) call cq_abort(' ERROR in recv_array_BtoG  : mx_recv_call,isend =',&
                  comm_naba_blocks_of_atoms%mx_recv_call,isend)
          call make_table(mynode,nnd_rem,npair,isend,iprim)
          call start_timer(tmr_std_allocation)
          deallocate(irecv_array,STAT=stat)
          if(stat/=0) call cq_abort("Error deallocating irecv_array in recv_array_BtoG: ",recv_size)
          call stop_timer(tmr_std_allocation)

       enddo
    endif

    return
  end subroutine recv_array_BtoG

  !-----------------------------------------------------------------------
  !sbrt make_table
  !  MAIN PART of this module
  !
  !  Modifications:
  !   2016/07/20 16:30 nakata
  !    Renamed naba_atm -> naba_atoms_of_blocks
  !   2016/08/01 17:30 nakata
  !    Introduced atomf instead of sf
  !-----------------------------------------------------------------------
  subroutine make_table(mynode,nnd_rem,npair,isend,iprim)

    use datatypes
    use numbers,       ONLY:very_small
    use global_module, ONLY:rcellx,rcelly,rcellz,x_atom_cell,y_atom_cell,z_atom_cell
    use group_module,  ONLY:parts,blocks
    use primary_module,ONLY:domain,bundle
    use cover_module,  ONLY:DCS_parts
    use global_module, ONLY: id_glob
    use atoms, ONLY: atoms_on_node, atom_number_on_node, node_doing_atom
    use GenComms, ONLY: cq_abort

    implicit none

    !Dummy Arguments
    integer,intent(in)::mynode,nnd_rem,npair,isend,iprim
    integer ::ncoverz,ncoveryz
    ! 28/Jul/2000
    integer :: ncoverx_add,ncovery_add,ncoverz_add
    real(double) :: dcellx_part,dcelly_part,dcellz_part
    real(double) :: dcellx_block,dcelly_block,dcellz_block
    !integer :: ipair,ind_block,ifind_block,iblock,ind_ib
    integer :: ipair,ind_block,ifind_block
    real(double) :: xmin,ymin,zmin
    integer :: ifind_part,ipart,jpart,ind_qart,nxp,nyp,nzp
    real(double) :: xmin_p,ymin_p,zmin_p
    integer :: nx,ny,nz,offset
    integer :: ind_part,ibegin,iend,ifind_pair,ii,j,ia1,ia2,iprim2
    real(double) :: dist,xatom,yatom,zatom,xx,yy,zz,dist_min
    integer :: ix,iy,iz
    integer :: irc,ierr

    ncoverx_add=DCS_parts%ncoverx
    ncovery_add=DCS_parts%ncovery
    ncoverz_add=DCS_parts%ncoverz

    ncoverz=2*DCS_parts%ncoverz+1
    ncoveryz=(2*DCS_parts%ncovery+1)*ncoverz

    dcellx_part=rcellx/parts%ngcellx ; dcellx_block=rcellx/blocks%ngcellx
    dcelly_part=rcelly/parts%ngcelly ; dcelly_block=rcelly/blocks%ngcelly
    dcellz_part=rcellz/parts%ngcellz ; dcellz_block=rcellz/blocks%ngcellz

    if(npair < 1) call cq_abort('npair in make_table <1 ?? : npair = ',nnd_rem,iprim)
    do ipair=1,npair
       ind_block=naba_blk_rem(ipair)
       ifind_block=0
       ifind_block=blocks%i_cc2seq(ind_block)
       if(ifind_block <= 0 .or. ifind_block > domain%groups_on_node) &
            call cq_abort('No block found in make_table: ',ipair,ifind_block)
       xmin=(domain%idisp_primx(ifind_block)+domain%nx_origin-1)*dcellx_block
       ymin=(domain%idisp_primy(ifind_block)+domain%ny_origin-1)*dcelly_block
       zmin=(domain%idisp_primz(ifind_block)+domain%nz_origin-1)*dcellz_block 
       ifind_part=0 

       do ipart=1,naba_atoms_of_blocks(atomf)%no_of_part(ifind_block) !Loop over naba parts
          jpart=naba_atoms_of_blocks(atomf)%list_part(ipart,ifind_block) !NOPG in a cover
          ind_qart=DCS_parts%lab_cover(jpart)-1             !CC in a cover
          nxp= ind_qart/(DCS_parts%ncovery*DCS_parts%ncoverz)
          nyp= (ind_qart-nxp*DCS_parts%ncovery*DCS_parts%ncoverz)/DCS_parts%ncoverz
          nzp= ind_qart-nxp*DCS_parts%ncovery*DCS_parts%ncoverz-nyp*DCS_parts%ncoverz
          xmin_p= (nxp+DCS_parts%nx_origin-DCS_parts%nspanlx-1)*dcellx_part
          ymin_p= (nyp+DCS_parts%ny_origin-DCS_parts%nspanly-1)*dcelly_part
          zmin_p= (nzp+DCS_parts%nz_origin-DCS_parts%nspanlz-1)*dcellz_part

          ! We must check carefully whether the following scheme 
          ! to calculate the offsets is consistent with the one used in
          ! sending nodes (bundle responsible nodes) in subroutine 
          ! 'get_naba_BCSblk'.
          !  especially which (part or block) coordinates will be added by very_small.

          nx=anint((xmin-xmin_p+very_small)/dcellx_part)
          ny=anint((ymin-ymin_p+very_small)/dcelly_part)
          nz=anint((zmin-zmin_p+very_small)/dcellz_part)

          nx=nx+DCS_parts%ncoverx
          ny=ny+DCS_parts%ncovery
          nz=nz+DCS_parts%ncoverz

          offset=(nx-1)*ncoveryz+(ny-1)*ncoverz+nz
          if(offset == offset_rem(ipair)) then
             ifind_part=ipart
             exit
          endif
       enddo ! enddo loop of neighbour partitions

       if(ifind_part == 0) call cq_abort('Pair not found in make_table: ',ipair,nnd_rem)
       jpart=naba_atoms_of_blocks(atomf)%list_part(ifind_part,ifind_block)
       ind_part=DCS_parts%lab_cell(jpart)

       ibegin=naba_atoms_of_blocks(atomf)%ibegin_part(ifind_part,ifind_block)
       iend=ibegin+naba_atoms_of_blocks(atomf)%no_atom_on_part(ifind_part,ifind_block)-1
       ifind_pair=0

       do ii=ibegin,iend    ! find the pair (OLD: check global no. of atom)
          !               (NOW: check primary no. of atom)
          j=naba_atoms_of_blocks(atomf)%list_atom(ii,ifind_block) !part seq. no. of atom
          ia1=id_glob(parts%icell_beg(ind_part)+j-1) 
          iprim2=atom_number_on_node(ia1)
          !ia1: glob no. of j-th atm in (ind_part)
          ia2=atoms_on_node(iprim,nnd_rem)           

          if(iprim /= atom_number_on_node(ia2)) &
               write(io_lun,*) ' ERROR: IPRIM in make_table! '
          !write(io_lun,102) ii,j,ia1,ia2
          !102 format(6x,'ii,j,ia1,ia2 = ',2i4,2i6)
          ! glob no. of iprim-th atm on nnd_rem
          !if(ia1 == ia2) then
          !iprim2=bundle%nm_nodbeg(parts%i_cc2seq(ind_part))+j WRONG!!
          if(iprim == iprim2) then
             ifind_pair=ii
             exit
          endif
       enddo  ! end loop of finding the pair

       if(ifind_pair == 0) call cq_abort('Atom not found in make_table: ',ipair,nnd_rem)
       !table_BtoG_npair(isend)        =npair    !no of pair  (not necessary)
       ! isend: seq. no. of receiving <-iprim&comm_naba_blocks_of_atoms%no_send_node(iprim)
       comm_naba_blocks_of_atoms%table_blk(ipair,isend) =ifind_block !shows prim block    
       comm_naba_blocks_of_atoms%table_pair(ipair,isend)=ifind_pair  !shows pair(naba atom)
    enddo ! ipair in my domain
    return
  end subroutine make_table
  !end subroutine make_table_BtoG

!!****f* set_blipgrid_module/free_blipgrid *
!!
!!  NAME 
!!   free_blipgrid
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
!!   08:06, 08/01/2003 dave
!!  MODIFICATION HISTORY
!!   11:53, 04/02/2003 drb 
!!    Bug fix - removed extraneous bracket
!!   2008/05/16 ast
!!    Added timer
!!   2016/07/06 17:30 nakata
!!    Renamed comBG to comm_naba_blocks_of_atoms
!!   2016/07/14 16:30 nakata
!!    Renamed naba_blk_supp -> naba_blocks_of_atoms
!!   2016/07/20 16:30 nakata
!!    Renamed naba_atm -> naba_atoms_of_blocks, halo_atm -> halo_atoms_of_blocks
!!   2016/08/01 17:30 nakata
!!    Introduced atomf instead of sf and paof
!!  SOURCE
!!
  subroutine free_blipgrid

    use naba_blk_module, ONLY: dealloc_halo_atm, dealloc_naba_atm, dealloc_naba_blk, dealloc_comm_in_BtoG

!    call start_timer(tmr_std_allocation)
    call dealloc_halo_atm(halo_atoms_of_blocks(dens))
!    call dealloc_halo_atm(halo_atoms_of_blocks(paof))   !!! nakata2 --- delete this line later
    call dealloc_halo_atm(halo_atoms_of_blocks(nlpf))
!    call dealloc_halo_atm(halo_atoms_of_blocks(sf))   !!! nakata2 --- delete this line later
    call dealloc_halo_atm(halo_atoms_of_blocks(atomf))   !!! nakata2
    call dealloc_naba_atm(naba_atoms_of_blocks(dens))
!    call dealloc_naba_atm(naba_atoms_of_blocks(paof))   !!! nakata2 --- delete this line later
    call dealloc_naba_atm(naba_atoms_of_blocks(nlpf))
!    call dealloc_naba_atm(naba_atoms_of_blocks(sf))   !!! nakata2 --- delete this line later
    call dealloc_naba_atm(naba_atoms_of_blocks(atomf))   !!! nakata2
    call dealloc_naba_blk(naba_blocks_of_atoms)
    call dealloc_comm_in_BtoG(comm_naba_blocks_of_atoms)
!    call stop_timer(tmr_std_allocation)
    return
  end subroutine free_blipgrid
!!***
end module set_blipgrid_module
