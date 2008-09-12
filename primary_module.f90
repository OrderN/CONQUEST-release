! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module primary_module
! ------------------------------------------------------------------------------
! Code area 9: general
! ------------------------------------------------------------------------------

!!****h* Conquest/primary_module
!!  NAME
!!   primary_module
!!  PURPOSE
!!   Build primary sets - for details of primary sets, small groups etc see the notes in matmult.tex 
!!   and the TechnicalMatMult notes in the CQDocs module
!!  USES
!!   GenComms, basic_types, datatypes, global_module
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   24/11/99
!!  MODIFICATION HISTORY
!!   30/03/00 by D.R.Bowler 
!!    Update global variable locations
!!   19/04/00 by DRB
!!    Generalised the idea of a primary set (so that we can use it for
!!    grid points in a domain OR atoms in a bundle), created a derived 
!!    type and changed make_prim to use them.
!!   20/04/00 by DRB
!!    Corrected assignment of nm_nodbeg and nm_nodgroup to primary_set
!!    Removed type declaration and added need for basic_types module
!!   22/05/2001 dave
!!    Added ROBODoc header and species attributes
!!   20/06/2001 dave
!!    Added RCS Id and Log tags and used cq_abort
!!   31/05/2002 dave
!!    Bug fix for make_prim (from TM), added RCSid and more comments
!!   2008/02/06 08:32 dave
!!    Changed for output to file not stdout
!!  SOURCE
!!
module primary_module

  ! Module usage 
  use datatypes
  use global_module, ONLY: io_lun
  use basic_types
  use GenComms, ONLY: cq_abort

  implicit none
  save

  ! Used for identifying object files (RCS ident command)
  character(len=80) :: RCSid = "$Id$"

  type(primary_set) :: domain ! Integration grid points
  type(primary_set) :: bundle ! Atoms
!!***

contains

!!****f* primary_module/make_prim *
!!
!!  NAME 
!!   make_prim
!!  USAGE
!!   make_prim(groups, prim[, m_id_glob, x_mem_cell, y_mem_cell, z_mem_cell, spec])
!!   make_prim(group type, primary set type[, global id, positions and species of members of group])
!!  PURPOSE
!!   Find the origin, left span, group offsets, group contents and
!! (optionally) group positions and species of a primary set
!!  INPUTS
!!   type(group_set) :: groups - type of group the set is made from
!!   type(primary_set) :: prim - the primary set we are building
!!   integer, dimension(:), OPTIONAL :: m_id_glob - global numbering of members
!!   real(double), dimension(:), OPTIONAL :: x_mem_cell - positions of members
!!   real(double), dimension(:), OPTIONAL :: y_mem_cell - positions of members
!!   real(double), dimension(:), OPTIONAL :: z_mem_cell - positions of members
!!   integer, dimension(:), OPTIONAL :: spec - species of members
!!  USES
!!   basic_types, global_module
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   Sometime in 2000
!!  MODIFICATION HISTORY
!!   22/05/2001 dave
!!    Added ROBODoc header and species information
!!   20/06/2001 dave
!!    Removed MPI_Abort
!!   31/05/2002 dave
!!    Added important bug fix from Tsuyoshi (zero iprim_seq only if members)
!!  SOURCE
!!
  subroutine make_prim(groups,prim,myid,m_id_glob, &
       x_mem_cell,y_mem_cell,z_mem_cell,spec)

    ! Module usage  
    use global_module
    use basic_types

    implicit none

    ! Passed variables
    type(group_set), target :: groups
    type(primary_set) :: prim
    integer, intent(in) :: myid
    integer, dimension(:), intent(IN), OPTIONAL :: m_id_glob
    real(double), dimension(:), intent(IN), OPTIONAL :: x_mem_cell
    real(double), dimension(:), intent(IN), OPTIONAL :: y_mem_cell
    real(double), dimension(:), intent(IN), OPTIONAL :: z_mem_cell
    integer, dimension(:), intent(IN), OPTIONAL :: spec

    ! Local variables
    ! iprojx(ip),y,z: shadows of primary-set groups on cell axes
    integer :: iprojx(groups%mx_gedge)
    integer :: iprojy(groups%mx_gedge)
    integer :: iprojz(groups%mx_gedge)
    integer :: irc,ierr,ind_group0,nx,ny,nz,ng,ind_group,nx1,ny1,nz1,ni
    integer :: nnd
    real(double) :: dcellx,dcelly,dcellz
    real(double) :: xadd,yadd,zadd
    logical :: members

    ! First determine if we're building details of members of set
    members = (PRESENT(m_id_glob).AND.PRESENT(x_mem_cell).AND.&
         PRESENT(y_mem_cell).AND.PRESENT(z_mem_cell))
    nnd = myid+1
    prim%groups => groups
    prim%groups_on_node = groups%ng_on_node(nnd)
    ! --- determine origin group of primary cell ----------------------
    ind_group0 = groups%ngnode(groups%inode_beg(nnd))
    prim%nx_origin = 1+(ind_group0-1)/(groups%ngcelly*groups%ngcellz)
    prim%ny_origin = 1+(ind_group0-1- &
         (prim%nx_origin-1)*groups%ngcelly*groups%ngcellz)/groups%ngcellz
    prim%nz_origin = ind_group0- &
         (prim%nx_origin-1)*groups%ngcelly*groups%ngcellz- &
         (prim%ny_origin-1)*groups%ngcellz
    if(iprint_gen>2.AND.myid==0) write(io_lun,1) myid,prim%nx_origin,prim%ny_origin,prim%nz_origin
1   format(2x,'On processor ',i4,' the primary set origin is ',3i5)
    ! --- determine widths and left spans of primary shell ----------------
    iprojx=0
    iprojy=0
    iprojz=0
    if(groups%inode_beg(nnd)+groups%ng_on_node(nnd)-1>groups%mx_gcell) then
       call cq_abort('make_prim: too many groups ',&
            groups%inode_beg(nnd)+groups%ng_on_node(nnd)-1)
    endif
    do ng=1,groups%ng_on_node(nnd)
       ind_group=groups%ngnode(groups%inode_beg(nnd)+ng-1)
       nx=1+(ind_group-1)/(groups%ngcelly*groups%ngcellz)
       ny=1+(ind_group-1-(nx-1)*groups%ngcelly*groups%ngcellz)/groups%ngcellz
       nz=ind_group-(nx-1)*groups%ngcelly*groups%ngcellz-(ny-1)*groups%ngcellz
       iprojx(nx)=1
       iprojy(ny)=1
       iprojz(nz)=1
    enddo
    call calliper(groups%mx_gedge,groups%ngcellx,prim%nx_origin, &
         iprojx,prim%nw_primx,prim%nleftx)
    call calliper(groups%mx_gedge,groups%ngcelly,prim%ny_origin, &
         iprojy,prim%nw_primy,prim%nlefty)
    call calliper(groups%mx_gedge,groups%ngcellz,prim%nz_origin, &
         iprojz,prim%nw_primz,prim%nleftz)
    ! --- determine offsets of groups from primary-set origin ---------
    do ng=1,prim%groups_on_node
       ind_group=groups%ngnode(groups%inode_beg(nnd)+ng-1)
       nx=1+(ind_group-1)/(groups%ngcelly*groups%ngcellz)
       ny=1+(ind_group-1-(nx-1)*groups%ngcelly*groups%ngcellz)/groups%ngcellz
       nz=ind_group-(nx-1)*groups%ngcelly*groups%ngcellz-(ny-1)*groups%ngcellz
       prim%idisp_primx(ng)=mod(nx+prim%nleftx-prim%nx_origin+groups%ngcellx, &
            groups%ngcellx)-prim%nleftx
       prim%idisp_primy(ng)=mod(ny+prim%nlefty-prim%ny_origin+groups%ngcelly, &
            groups%ngcelly)-prim%nlefty
       prim%idisp_primz(ng)=mod(nz+prim%nleftz-prim%nz_origin+groups%ngcellz, &
            groups%ngcellz)-prim%nleftz
    enddo
    ! --- analyse atoms numbers and positions in primary cell -------------
    dcellx=rcellx/groups%ngcellx
    dcelly=rcelly/groups%ngcelly
    dcellz=rcellz/groups%ngcellz
    if(members) prim%iprim_seq = 0
    prim%n_prim=0
    prim%nm_nodbeg(1)=1
    do ng=1,groups%ng_on_node(nnd)
       ind_group=groups%ngnode(groups%inode_beg(nnd)+ng-1)
       nx=1+(ind_group-1)/(groups%ngcelly*groups%ngcellz)
       ny=1+(ind_group-1-(nx-1)*groups%ngcelly*groups%ngcellz)/groups%ngcellz
       nz=ind_group-(nx-1)*groups%ngcelly*groups%ngcellz-(ny-1)*groups%ngcellz
       nx1=prim%nx_origin+prim%idisp_primx(ng)
       ny1=prim%ny_origin+prim%idisp_primy(ng)
       nz1=prim%nz_origin+prim%idisp_primz(ng)
       xadd=real(nx1-nx,double)*dcellx
       yadd=real(ny1-ny,double)*dcelly
       zadd=real(nz1-nz,double)*dcellz
       prim%nm_nodgroup(ng)=groups%nm_group(ind_group)
       if(prim%nm_nodgroup(ng).gt.0) then
          do ni=1,prim%nm_nodgroup(ng)
             prim%n_prim=prim%n_prim+1
             if(prim%n_prim.gt.prim%mx_iprim) then
                call cq_abort('make_prim: too many atoms ', &
                     prim%n_prim+prim%nm_nodgroup(ng)-ni,prim%mx_iprim)
             endif
             if(members) then
                prim%iprim_seq(prim%n_prim) = ni
                prim%ig_prim(prim%n_prim)= &
                     m_id_glob(groups%icell_beg(ind_group)+ni-1)
                prim%xprim(prim%n_prim)= &
                     x_mem_cell(groups%icell_beg(ind_group)+ni-1)+xadd
                prim%yprim(prim%n_prim)= &
                     y_mem_cell(groups%icell_beg(ind_group)+ni-1)+yadd
                prim%zprim(prim%n_prim)= &
                     z_mem_cell(groups%icell_beg(ind_group)+ni-1)+zadd
                prim%species(prim%n_prim)= &
                     spec(groups%icell_beg(ind_group)+ni-1)
                if(iprint_gen>4.AND.myid==0) write(io_lun,fmt='(2x,"Prim atom: ",i4," position: ",3f8.3)') prim%n_prim, &
                     prim%xprim(prim%n_prim),prim%yprim(prim%n_prim),prim%zprim(prim%n_prim)
             endif
          enddo
       endif
       if(ng.lt.groups%ng_on_node(nnd)) then
          prim%nm_nodbeg(ng+1)=prim%nm_nodbeg(ng)+prim%nm_nodgroup(ng)
       endif
    enddo
    return
  end subroutine make_prim
!!***

!!****f* primary_module/calliper *
!!
!!  NAME 
!!   calliper
!!  USAGE
!! 
!!  PURPOSE
!!   Determines width and left-hand span of primary
!!   set in one of the cartesian directions.
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   M.J.Gillan/D.R.Bowler
!!  CREATION DATE
!! 
!!  MODIFICATION HISTORY
!!   22/05/2001 dave
!!    Added ROBODoc header
!!  SOURCE
!!
  subroutine calliper(mx_gedge,ngcell,n_origin,iproj,nw_prim,nleft)

    implicit none

    ! Passed variables
    integer :: mx_gedge,ngcell,n_origin,nw_prim,nleft
    integer :: iproj(mx_gedge)

    ! Local variables
    integer :: irc,ierr,n_in_gap,ngap_beg,i,im,n_in_gap1,i1,i_left

    ! Find the gap
    n_in_gap=0
    ngap_beg=1
    do i=1,ngcell
       im=1+mod(i+ngcell-2,ngcell)
       if((iproj(im).eq.1).and.(iproj(i).eq.0)) then
          n_in_gap1=0
          i1=i
          do while(iproj(i1).eq.0)
             n_in_gap1=n_in_gap1+1
             i1=1+mod(i1,ngcell)
          enddo
          if(n_in_gap1.gt.n_in_gap) then
             n_in_gap=n_in_gap1
             ngap_beg=i
          endif
       endif
    enddo
    nw_prim=ngcell-n_in_gap
    i_left=1+mod(ngap_beg+n_in_gap-1,ngcell)
    nleft=mod(n_origin-i_left+ngcell,ngcell)
    return
  end subroutine calliper
!!***

!!****f* primary_module/allocate_primary_set *
!!
!!  NAME 
!!   allocate_primary_set
!!  USAGE
!! 
!!  PURPOSE
!!   Allocates memory to a primary set structure
!!   N.B.: the mx_iprim variable is OPTIONAL - typically only used for
!!   the atomic positions when groups are partitions
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
!!   22/05/2001 dave
!!    ROBODoc header
!!   20/06/2001 dave
!!    Added cq_abort
!!  SOURCE
!!
  subroutine allocate_primary_set(prim,members)

    ! Module usage
    use datatypes
    use basic_types

    ! Passed variables
    type(primary_set) :: prim
    logical :: members

    ! Local variables
    integer :: stat

    allocate(prim%idisp_primx(prim%mx_ngonn),STAT=stat)
    if(stat/=0) then
       call cq_abort('alloc_prim: error(1) idisp_primx')
    endif
    allocate(prim%idisp_primy(prim%mx_ngonn),STAT=stat)
    if(stat/=0) then
       call cq_abort('alloc_prim: error(2) idisp_primy')
    endif
    allocate(prim%idisp_primz(prim%mx_ngonn),STAT=stat)
    if(stat/=0) then
       call cq_abort('alloc_prim: error(3) idisp_primz')
    endif
    allocate(prim%nm_nodgroup(prim%mx_ngonn),STAT=stat)
    if(stat/=0) then
       call cq_abort('alloc_prim: error(4) nm_nodgroup')
    endif
    allocate(prim%nm_nodbeg(prim%mx_ngonn),STAT=stat)
    if(stat/=0) then
       call cq_abort('alloc_prim: error(5) nm_nodbeg')
    endif
    if(members) then
       allocate(prim%iprim_seq(prim%mx_iprim),STAT=stat)
       if(stat/=0) then
          call cq_abort('alloc_prim: error(6) iprim_seq')
       endif
       allocate(prim%ig_prim(prim%mx_iprim),STAT=stat)
       if(stat/=0) then
          call cq_abort('alloc_prim: error(7) ig_prim')
       endif
       allocate(prim%xprim(prim%mx_iprim),STAT=stat)
       if(stat/=0) then
          call cq_abort('alloc_prim: error(8) xprim')
       endif
       allocate(prim%yprim(prim%mx_iprim),STAT=stat)
       if(stat/=0) then
          call cq_abort('alloc_prim: error(9) yprim')
       endif
       allocate(prim%zprim(prim%mx_iprim),STAT=stat) 
       if(stat/=0) then
          call cq_abort('alloc_prim: error(10) zprim')
       endif
       allocate(prim%species(prim%mx_iprim),STAT=stat) 
       if(stat/=0) then
          call cq_abort('alloc_prim: error(11) species')
       endif
    else
       nullify(prim%iprim_seq,prim%ig_prim, &
            prim%xprim,prim%yprim,prim%zprim,prim%species)
    endif
    return
  end subroutine allocate_primary_set
!!***

!!****f* primary_module/deallocate_primary_set *
!!
!!  NAME 
!!   deallocate_primary_set
!!  USAGE
!! 
!!  PURPOSE
!!   Deallocates the primary set structure
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
!!   22/05/2001 dave
!!    ROBODoc header
!!   20/06/2001 dave
!!    Added cq_abort
!!  SOURCE
!!
  subroutine deallocate_primary_set(prim)

    ! Module usage
    use datatypes
    use basic_types

    implicit none

    ! Passed variables
    type(primary_set) :: prim

    ! Local variables
    integer :: stat

    deallocate(prim%species,prim%zprim,prim%yprim,prim%xprim, &
         prim%ig_prim,prim%iprim_seq, &
         prim%nm_nodbeg,prim%nm_nodgroup, &
         prim%idisp_primz,prim%idisp_primy,prim%idisp_primx,STAT=stat)
    if(stat/=0) then
       call cq_abort('dealloc_prim: error deallocating')
    endif
    return
  end subroutine deallocate_primary_set
!!***
end module primary_module
