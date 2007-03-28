! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id: cover_module.f90,v 1.6 2003/10/22 15:06:25 mjg Exp $
! ------------------------------------------------------------------------------
! Module cover_module
! ------------------------------------------------------------------------------
! Code area 8: indexing
! -------------------------------------------------------------------------------

!!****h* Conquest/cover_module
!!  NAME
!!   cover_module
!!  PURPOSE
!!   This module holds variables and routines related to covering
!!   sets.  A covering set is (at the moment) the smallest orthorhombic
!!   collection of small groups that will fully enclose a primary set
!!   out to a specified radius.  It's used for searching for neighbours
!!   and communicating details about small groups between processors.
!!
!!   More details (exhaustive details !) can be found in the Conquest
!!   notes "Matrix multiplication in Conquest: Practical details" 
!!   which are in the Conquest documentation repository (stored
!!   in the directory TechnicalMatMult/)
!!  USES
!!   GenComms, basic_types, datatypes, global_module
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   24/11/99
!!  MODIFICATION HISTORY
!!   19/02/00 by D.R.Bowler
!!    Implementing changes for transposes (adding ncover_rem and
!!    routines to find iprim_part)
!!   03/03/00 by DRB
!!    Finishing up transposes
!!   29/03/00 by DRB 
!!    Creating a generalised cover type and altering make_gcs to use it
!!   11/04/00 by DRB
!!    Correcting cover type
!!   19/04/00 by DRB
!!    Adding general group and primary set types to make GCS totally 
!!    general
!!   20/04/00 by DRB
!!    Removed type declaration and placed in basic_types module
!!    Continued generalisation of make_cs
!!    Removed GOTOs from indexx (replaced with do...while)
!!   21/06/2001 dave
!!    Added ROBODoc headers and removed stop statements and MPI_Abort
!!    statements throughout
!!   19/04/2002 drb
!!    Corrected reference to notes and added USE field for whole module
!!   17/06/2002 dave
!!    Added RCS Id tag and improved headers a little
!!  SOURCE
!!
module cover_module

  ! Module use
  use datatypes
  use basic_types

  implicit none
  save

  ! DCS = domain CS; BCS = bundle CS
  type(cover_set) :: DCS_blocks
  type(cover_set) :: DCS_parts
  type(cover_set) :: BCS_blocks
  type(cover_set) :: BCS_parts
  type(cover_set) :: ewald_CS

  ! -------------------------------------------------------
  ! RCS ident string for object file id
  ! -------------------------------------------------------
  character(len=80), private :: RCSid = "$Id: cover_module.f90,v 1.6 2003/10/22 15:06:25 mjg Exp $"

!!***
contains

!!****f* cover_module/make_cs *
!!
!!  NAME 
!!   make_cs - makes a covering set
!!  USAGE
!! 
!!  PURPOSE
!!   sbrt make_cs: create a covering set out of groups given by the 
!!   group_set variable, for a cutoff gr_rcut around the primary set 
!!   defined by the widths (nw_x,y,z) and left-spans (nlx,y,z).
!!   N.B. The widths and left spans must have been converted into units
!!   of the CS groups !
!!   The origin is also passed; this is a problem is we're making a
!!   CS for one group type out of the other.  It MUST be the group 
!!   nearest the origin partition of the primary set.
!!   The results are put into the derived type set.
!!   The logical variable members decides whether or not the 
!!   x,y,zcover variables are created (NOT for the cross-group CS)
!!
!!   The notes above are a little out of date - the widths and left-spans
!!   are now part of the primary set type, and conversion is done by 
!!   convert_primary (below)
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   24/11/99
!!  MODIFICATION HISTORY
!!   21/06/2001 dave
!!    Added ROBODoc header and cq_abort
!!  SOURCE
!!
  subroutine make_cs(myid,gr_rcut,set,groups,prim, mx_mcell,x_mem_cell,y_mem_cell,z_mem_cell)

    ! Module usage
    use global_module
    use basic_types
    use GenComms, ONLY: cq_abort, inode, ionode
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_int, type_dbl

    implicit none

    ! Passed variables
    type(cover_set) :: set
    type(group_set), target   :: groups
    type(primary_set), target :: prim
    integer, intent(in) :: myid
    integer, OPTIONAL :: mx_mcell
    real(double), dimension(:), intent(IN), OPTIONAL :: x_mem_cell
    real(double), dimension(:), intent(IN), OPTIONAL :: y_mem_cell
    real(double), dimension(:), intent(IN), OPTIONAL :: z_mem_cell
    real(double) :: gr_rcut ! Cutoff radius for CS

    ! Local variables
    ! Number of periodic images of a group
    integer :: nrepx(groups%mx_gedge)
    integer :: nrepy(groups%mx_gedge)
    integer :: nrepz(groups%mx_gedge)
    ! x,y,z numbering of CS groups (for CC labels)
    integer, allocatable, dimension(:)  :: nx_in_cover
    integer, allocatable, dimension(:)  :: ny_in_cover
    integer, allocatable, dimension(:)  :: nz_in_cover
    ! Variables for irreducible CS
    integer :: ind_min(groups%mx_gcell)
    integer :: ngcx_min(groups%mx_gcell)
    integer :: ngcy_min(groups%mx_gcell)
    integer :: ngcz_min(groups%mx_gcell)
    integer :: min_sort(groups%mx_gcell)
    ! Converted units of primary cell size
    integer :: nx_o,ny_o,nz_o
    ! General scalars
    integer :: ng_in_cell,nmodx,nmody,nmodz
    integer :: noccx,nremx,minx,ngcx,noccy,nremy,miny,ngcy,noccz,nremz
    integer :: minz,ngcz,ng_in_min,ind,nqx,nqy,nqz,ind_qart,ino,ind_cover
    integer :: nrx,nry,nrz,nsx,nsy,nsz,ni,nnd,irc,ierr,stat,pr,i,j
    integer :: ind_part,np, nm_in_cover
    real(double) :: dcellx,dcelly,dcellz,xadd,yadd,zadd
    logical :: members

    ! This determines whether or not we find information about all 
    ! the members of the groups in the covering set
    members = (PRESENT(mx_mcell).AND.PRESENT(x_mem_cell).AND.&
         PRESENT(y_mem_cell).AND.PRESENT(z_mem_cell))  
    nnd = myid+1
    ! Start the derived type - assign useful variables
    set%rcut = gr_rcut     ! Cutoff around primary set
    set%groups => groups   ! Groups out of which to build CS
    set%prim => prim       ! Primary set around which to build CS
    ! Find the origin, width and left span of the CS
    call convert_primary(groups,prim,set,nx_o,ny_o,nz_o)
    ! This gives the total number of groups
    set%ng_cover=set%ncoverx*set%ncovery*set%ncoverz
    set%mx_gcover = set%ng_cover
    call allocate_cs(set,numprocs,members)
    allocate(nx_in_cover(set%ng_cover),ny_in_cover(set%ng_cover),nz_in_cover(set%ng_cover),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating nx_in_cover: ",set%ng_cover,stat)
    call reg_alloc_mem(area_index, 3*set%ng_cover,type_int)
    !if(set%ng_cover>set%mx_gcover) then
    !   call cq_abort('make_cs: too many groups in GCS',&
    !        set%ng_cover,set%mx_gcover)
    !endif
    ! Conversion factors from unit cell lengths->groups
    dcellx=rcellx/real(groups%ngcellx,double)
    dcelly=rcelly/real(groups%ngcelly,double)
    dcellz=rcellz/real(groups%ngcellz,double)
    ! Fully explained in notes mentioned above
    nmodx=((groups%ngcellx+set%nspanlx-1)/groups%ngcellx)*groups%ngcellx
    nmody=((groups%ngcelly+set%nspanly-1)/groups%ngcelly)*groups%ngcelly
    nmodz=((groups%ngcellz+set%nspanlz-1)/groups%ngcellz)*groups%ngcellz
    ! create apparatus for periodic boundary conditions ---------------
    ! ... x-direction
    noccx=set%ncoverx/groups%ngcellx
    nremx=set%ncoverx-noccx*groups%ngcellx
    minx=min(set%ncoverx,groups%ngcellx)
    if(minx>groups%mx_gedge) then
       call cq_abort('make_cs: too many groups in x-edge',minx)
    endif
    do ngcx=1,minx
       if(ngcx<=nremx) then
          nrepx(ngcx)=noccx+1
       else
          nrepx(ngcx)=noccx
       endif
    enddo
    ! ... y-direction
    noccy=set%ncovery/groups%ngcelly
    nremy=set%ncovery-noccy*groups%ngcelly
    miny=min(set%ncovery,groups%ngcelly)
    if(miny>groups%mx_gedge) then
       call cq_abort('make_cs: too many groups in y-edge',miny)
    endif
    do ngcy=1,miny
       if(ngcy<=nremy) then
          nrepy(ngcy)=noccy+1
       else
          nrepy(ngcy)=noccy
       endif
    enddo
    ! ... z-direction
    noccz=set%ncoverz/groups%ngcellz
    nremz=set%ncoverz-noccz*groups%ngcellz
    minz=min(set%ncoverz,groups%ngcellz)
    if(minz>groups%mx_gedge) then
       call cq_abort('make_cs: too many groups in z-edge',minz)
    endif
    do ngcz=1,minz
       if(ngcz<=nremz) then
          nrepz(ngcz)=noccz+1
       else
          nrepz(ngcz)=noccz
       endif
    enddo
    ! go over groups in GCS periodic-irreducible set, calculating
    ! simulation-cell (node-order, home-start) label of each 
    ng_in_cell = groups%ngcellx*groups%ngcelly*groups%ngcellz
    ng_in_min = minx*miny*minz
    ind=0
    do ngcx=1,minx
       do ngcy=1,miny
          do ngcz=1,minz
             ind=ind+1
             nqx=1+mod(nx_o+ngcx-set%nspanlx-2+nmodx,groups%ngcellx)
             nqy=1+mod(ny_o+ngcy-set%nspanly-2+nmody,groups%ngcelly)
             nqz=1+mod(nz_o+ngcz-set%nspanlz-2+nmodz,groups%ngcellz)
             ind_qart=(nqx-1)*groups%ngcelly*groups%ngcellz+(nqy-1)*groups%ngcellz+nqz
             ino=groups%inv_ngnode(ind_qart)
             ind_min(ind)=1+mod(ino-groups%inode_beg(nnd)+ ng_in_cell,ng_in_cell)
             ngcx_min(ind)=ngcx
             ngcy_min(ind)=ngcy
             ngcz_min(ind)=ngcz
          enddo
       enddo
    enddo
    ! sort minimum CS by nodes 
    call indexx(groups%mx_gcell,ng_in_min,ind_min,min_sort)
    ! go over all GCS groups in node-periodic-grouped order 
    ind_cover=0
    do ind=1,ng_in_min
       ngcx=ngcx_min(min_sort(ind))
       ngcy=ngcy_min(min_sort(ind))
       ngcz=ngcz_min(min_sort(ind))
       do nrx=1,nrepx(ngcx)
          do nry=1,nrepy(ngcy)
             do nrz=1,nrepz(ngcz)
                ind_cover=ind_cover+1
                nx_in_cover(ind_cover)=ngcx-1-set%nspanlx+(nrx-1)*groups%ngcellx
                ny_in_cover(ind_cover)=ngcy-1-set%nspanly+(nry-1)*groups%ngcelly
                nz_in_cover(ind_cover)=ngcz-1-set%nspanlz+(nrz-1)*groups%ngcellz
             enddo
          enddo
       enddo
    enddo
    ! make covering set for primary set 
    !if(members) set%icover_ibeg(1)=1
    if(members) nm_in_cover = 0
    do ind_cover=1,set%ng_cover
       nsx=nx_in_cover(ind_cover)
       nsy=ny_in_cover(ind_cover)
       nsz=nz_in_cover(ind_cover)
       nqx=1+mod(nx_o+nsx+nmodx-1,groups%ngcellx)
       nqy=1+mod(ny_o+nsy+nmody-1,groups%ngcelly)
       nqz=1+mod(nz_o+nsz+nmodz-1,groups%ngcellz)
       xadd=(nx_o+nsx-nqx)*dcellx
       yadd=(ny_o+nsy-nqy)*dcelly
       zadd=(nz_o+nsz-nqz)*dcellz
       ind_qart=(nqx-1)*groups%ngcelly*groups%ngcellz+&
            (nqy-1)*groups%ngcellz+nqz
       set%lab_cell(ind_cover)=ind_qart
       set%lab_cover(ind_cover)=&
            (nsx+set%nspanlx)*set%ncovery*set%ncoverz+&
            (nsy+set%nspanly)*set%ncoverz+(nsz+set%nspanlz)+1
       !!  TM  temporary ????
       set%inv_lab_cover(set%lab_cover(ind_cover))=ind_cover
       !!  TM  temporary ????
       if(members) nm_in_cover = nm_in_cover + groups%nm_group(ind_qart)
       !if(members.AND.groups%nm_group(ind_qart)>0) then
       !   set%n_ing_cover(ind_cover)=groups%nm_group(ind_qart)
       !   if(set%icover_ibeg(ind_cover)+groups%nm_group(ind_qart)-1> &
       !        set%mx_mcover) then
       !      call cq_abort('make_cs: xcover dimension exceeded',& 
       !           set%icover_ibeg(ind_cover)+groups%nm_group(ind_qart)-1, &
       !           set%mx_mcover)
       !   endif
       !   if(groups%icell_beg(ind_qart)+groups%nm_group(ind_qart)-1> &
       !        mx_mcell) then
       !      call cq_abort('make_cs: x_atom_cover dim. exceeded', &
       !           groups%icell_beg(ind_qart)+groups%nm_group(ind_qart)-1)
       !   endif
       !   do ni=1,groups%nm_group(ind_qart)
       !      set%xcover(set%icover_ibeg(ind_cover)+ni-1)= &
       !           x_mem_cell(groups%icell_beg(ind_qart)+ni-1)+xadd
       !      set%ycover(set%icover_ibeg(ind_cover)+ni-1)= &
       !           y_mem_cell(groups%icell_beg(ind_qart)+ni-1)+yadd
       !      set%zcover(set%icover_ibeg(ind_cover)+ni-1)= &
       !           z_mem_cell(groups%icell_beg(ind_qart)+ni-1)+zadd
       !   enddo
       !   if(ind_cover.lt.set%ng_cover) then
       !      set%icover_ibeg(ind_cover+1)=set%icover_ibeg(ind_cover)+ &
       !           groups%nm_group(ind_qart)
       !   endif
       !else if(members) then
       !   if(ind_cover < set%ng_cover) then   !TM 26/Jun/2003
       !      set%icover_ibeg(ind_cover+1)=set%icover_ibeg(ind_cover)+ &
       !           groups%nm_group(ind_qart)
       !   endif  ! (ind_cover < set%ng_cover)  !TM 26/Jun/2003
       !   set%n_ing_cover(ind_cover)=0
       !endif ! members.AND.nm_group>0
    enddo
    if(members) then ! Now generate member information
       if(inode==ionode.AND.iprint_index>1) write(*,*) 'Members in covering set: ',nm_in_cover
       if(members) set%icover_ibeg(1)=1
       set%mx_mcover = nm_in_cover
       allocate(set%xcover(set%mx_mcover),STAT=stat)
       if(stat/=0) call cq_abort('Error allocating cover set members x: ',set%mx_mcover,stat)
       allocate(set%ycover(set%mx_mcover),STAT=stat)
       if(stat/=0) call cq_abort('Error allocating cover set members y: ',set%mx_mcover,stat)
       allocate(set%zcover(set%mx_mcover),STAT=stat)
       if(stat/=0) call cq_abort('Error allocating cover set members z: ',set%mx_mcover,stat)
       call reg_alloc_mem(area_index, 3*set%mx_mcover,type_dbl)
       do ind_cover=1,set%ng_cover
          nsx=nx_in_cover(ind_cover)
          nsy=ny_in_cover(ind_cover)
          nsz=nz_in_cover(ind_cover)
          nqx=1+mod(nx_o+nsx+nmodx-1,groups%ngcellx)
          nqy=1+mod(ny_o+nsy+nmody-1,groups%ngcelly)
          nqz=1+mod(nz_o+nsz+nmodz-1,groups%ngcellz)
          xadd=(nx_o+nsx-nqx)*dcellx
          yadd=(ny_o+nsy-nqy)*dcelly
          zadd=(nz_o+nsz-nqz)*dcellz
          ind_qart= set%lab_cell(ind_cover)
          if(groups%nm_group(ind_qart)>0) then
             set%n_ing_cover(ind_cover)=groups%nm_group(ind_qart)
             if(groups%icell_beg(ind_qart)+groups%nm_group(ind_qart)-1> mx_mcell) then
                call cq_abort('make_cs: x_atom_cover dim. exceeded', &
                     groups%icell_beg(ind_qart)+groups%nm_group(ind_qart)-1)
             endif
             do ni=1,groups%nm_group(ind_qart)
                set%xcover(set%icover_ibeg(ind_cover)+ni-1)= &
                     x_mem_cell(groups%icell_beg(ind_qart)+ni-1)+xadd
                set%ycover(set%icover_ibeg(ind_cover)+ni-1)= &
                     y_mem_cell(groups%icell_beg(ind_qart)+ni-1)+yadd
                set%zcover(set%icover_ibeg(ind_cover)+ni-1)= &
                     z_mem_cell(groups%icell_beg(ind_qart)+ni-1)+zadd
             enddo
             if(ind_cover.lt.set%ng_cover) then
                set%icover_ibeg(ind_cover+1)=set%icover_ibeg(ind_cover)+ &
                     groups%nm_group(ind_qart)
             endif
          else 
             if(ind_cover < set%ng_cover) then   !TM 26/Jun/2003
                set%icover_ibeg(ind_cover+1)=set%icover_ibeg(ind_cover)+ &
                     groups%nm_group(ind_qart)
             endif  ! (ind_cover < set%ng_cover)  !TM 26/Jun/2003
             set%n_ing_cover(ind_cover)=0
          endif ! nm_group>0
       enddo
    end if
    !if(members) then ! Added by TM 17/07/00
    !   do np=1,set%ng_cover  ! Loop over partitions in GCS
    !      if(set%n_ing_cover(np).gt.0) then  ! Are there atoms ?
    !         if(set%icover_ibeg(np)+set%n_ing_cover(np)-1.gt.&
    !              set%mx_mcover) then
    !            call cq_abort('make_cs: member index error', &
    !                 set%icover_ibeg(np)+set%n_ing_cover(np)-1)
    !         endif
    !      endif
    !   enddo
    !endif ! (members)
    deallocate(nx_in_cover,ny_in_cover,nz_in_cover,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating nx_in_cover: ",set%ng_cover,stat)
    call reg_dealloc_mem(area_index, 3*set%ng_cover,type_int)
    return
  end subroutine make_cs
!!***

!!****f* cover_module/make_iprim *
!!
!!  NAME 
!!   make_iprim
!!  USAGE
!! 
!!  PURPOSE
!!   sbrt make_iprim: creates iprim_group, which gives the covering set CC 
!!   label for each primary set group - I'm pretty sure that this is 
!!   right, but it's a really hard problem
!!
!!   Definitely right - it was broken a while back, but now fixed 
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   24/11/99 ?
!!  MODIFICATION HISTORY
!!   21/06/2001 dave
!!    Added ROBODoc header and cq_abort
!!  SOURCE
!!
  subroutine make_iprim(set,prim,nnd)

    use basic_types
    use GenComms, ONLY: cq_abort
    use memory_module, ONLY: reg_alloc_mem, type_int
    use global_module, ONLY: area_index

    implicit none

    ! Passed variables
    type(cover_set) :: set
    type(primary_set) :: prim
    integer :: nnd

    ! Local variables
    integer :: pr,i,nsx,nsy,nsz,ind_group,j,irc,ierr, stat

    if(.NOT.ASSOCIATED(set%iprim_group)) then
       allocate(set%iprim_group(prim%mx_iprim),STAT=stat)
       if(stat/=0) call cq_abort('make_iprim: error allocating memory')
       call reg_alloc_mem(area_index,prim%mx_iprim,type_int)
    endif
    pr = 0 ! Indexes primary members
    do i=1,prim%groups_on_node
       nsx = prim%idisp_primx(i)
       nsy = prim%idisp_primy(i)
       nsz = prim%idisp_primz(i)
       ind_group = (nsx+set%nspanlx)*set%ncovery*set%ncoverz+&
            (nsy+set%nspanly)*set%ncoverz+(nsz+set%nspanlz)+1
       if(prim%nm_nodgroup(i) > 0) then
          do j=1,prim%nm_nodgroup(i)
             pr = pr+1
             if(pr>prim%mx_iprim) call cq_abort('make_iprim: error in pr index',pr)
             set%iprim_group(pr) = ind_group
          enddo
       endif ! (prim%nm_nodgroup(i) > 0) then
    enddo
  end subroutine make_iprim
!!***

!!****f* cover_module/send_ncover *
!!
!!  NAME 
!!   send_ncover
!!  USAGE
!! 
!!  PURPOSE
!!   Once make_cs has been run on every processor, this routine
!!   sends sizes of the CS to all other processors (and gathers theirs)
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   24/11/99 ?
!!  MODIFICATION HISTORY
!!   21/06/2001 dave
!!    ROBODoc header and changed to gsum
!!  SOURCE
!!
! =====================================================================
!   sbrt send_ncover: Called once at the start of a run (after make_gcs), 
!   to ensure that all nodes know all other nodes GCS sizes (ncoverx)
!   for transpose work
! ---------------------------------------------------------------------
  subroutine send_ncover(set,nnd)

    ! Module declaration
    use global_module, ONLY: numprocs
    use GenComms, ONLY: gsum, cq_abort

    implicit none

    ! Passed variables
    type(cover_set) :: set
    integer :: nnd ! Node no. starting from one

    ! Local variables
    integer, allocatable, dimension(:) :: ncovers
    integer :: ierr,i,stat

    allocate(ncovers(3*numprocs),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating ncovers: ",stat,numprocs)
    set%ncover_rem = 0
    set%ncover_rem(3*(nnd-1)+1) = set%ncoverx
    set%ncover_rem(3*(nnd-1)+2) = set%ncovery
    set%ncover_rem(3*(nnd-1)+3) = set%ncoverz
    call gsum(set%ncover_rem, 3*numprocs)
    deallocate(ncovers,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating ncovers: ",stat,numprocs)
    return
  end subroutine send_ncover
!!***

!!****f* cover_module/convert_primary *
!!
!!  NAME 
!!   convert_primary
!!  USAGE
!! 
!!  PURPOSE
!!   Given a primary set and a cutoff, creates the widths and left-spans
!!   of the CS in the specified groups and returns the origin of the CS
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   24/11/99 ?
!!  MODIFICATION HISTORY
!!   21/06/2001 dave
!!    Added ROBODoc header and used cq_abort
!!  SOURCE
!!
  subroutine convert_primary(groups,prim,set,ncx_o,ncy_o,ncz_o)

    use datatypes
    use global_module
    use basic_types
    use GenComms, ONLY: cq_abort

    implicit none

    ! Passed variables
    type(group_set)   :: groups
    type(primary_set) :: prim
    type(cover_set)   :: set
    integer :: ncx_o,ncy_o,ncz_o                ! Origin for cover set

    ! Local variables
    real(double), parameter :: eps = 0.00000001_double

    integer :: irc,ierr  
    real(double) :: ro_x,ro_y,ro_z              ! Origin in reals
    real(double) :: ro_cx,ro_cy,ro_cz           ! CS origin in reals
    real(double) :: dx,dy,dz,dcx,dcy,dcz        ! Sizes of groups in reals
    real(double) :: lhx,lhy,lhz,lhcx,lhcy,lhcz  ! Left hand corners
    real(double) :: rhx,rhy,rhz,rhcx,rhcy,rhcz  ! Right hand corners
    real(double) :: cslx,csly,cslz              ! Covering set left

    if(.NOT.associated(prim%groups)) then
       call cq_abort('convert_primary: primary set groups pointer not assoc')
    endif
    ! First, find sizes of groups in terms of unit cell
    dx = rcellx/dble(prim%groups%ngcellx)   ! Check for a more portable way !
    dy = rcelly/dble(prim%groups%ngcelly)
    dz = rcellz/dble(prim%groups%ngcellz)
    dcx = rcellx/dble(groups%ngcellx)
    dcy = rcelly/dble(groups%ngcelly)
    dcz = rcellz/dble(groups%ngcellz)
    !    write(*,*) 'dx, dcx: ',dx,dy,dz,dcx,dcy,dcz
    ! Convert origin of prim to reals (add eps to prevent ambiguity)
    ro_x = (dble(prim%nx_origin-1))*dx + eps 
    ro_y = (dble(prim%ny_origin-1))*dy + eps
    ro_z = (dble(prim%nz_origin-1))*dz + eps
    ! Convert this to number of CS groups
    ro_cx = ro_x/dcx
    ro_cy = ro_y/dcy
    ro_cz = ro_z/dcz
    ! Round and offset (partition 1,1,1 has LH corner at 0.0,0.0,0.0)
    ncx_o = 1+anint(ro_cx)  
    ncy_o = 1+anint(ro_cy)
    ncz_o = 1+anint(ro_cz)
    set%nx_origin = ncx_o
    set%ny_origin = ncy_o
    set%nz_origin = ncz_o
    ! Find CS origin in reals
    ro_cx = (dble(ncx_o-1))*dcx
    ro_cy = (dble(ncy_o-1))*dcy
    ro_cz = (dble(ncz_o-1))*dcz
    !    write(*,*) 'Origin: ',ro_x,ro_y,ro_z,ro_cx,ro_cy,ro_cz
    ! Start by finding the left and right hand corners of the CS
    lhx = dx*(prim%nx_origin-1-prim%nleftx)-set%rcut
    lhy = dy*(prim%ny_origin-1-prim%nlefty)-set%rcut
    lhz = dz*(prim%nz_origin-1-prim%nleftz)-set%rcut
    rhx = dx*(prim%nx_origin-1-prim%nleftx+prim%nw_primx)+set%rcut  
    rhy = dy*(prim%ny_origin-1-prim%nlefty+prim%nw_primy)+set%rcut
    rhz = dz*(prim%nz_origin-1-prim%nleftz+prim%nw_primz)+set%rcut
    ! The cover set left span must be offset from the CS origin
    cslx = (ro_cx - lhx)/dcx
    csly = (ro_cy - lhy)/dcy
    cslz = (ro_cz - lhz)/dcz
    ! Convert to integers
    set%nspanlx = ceiling(cslx)
    set%nspanly = ceiling(csly)
    set%nspanlz = ceiling(cslz)
    ! The width of the CS must be the distance from the LH corner to the 
    ! RH corner - this way we pick the SMALLEST CS needed.
    set%ncoverx = ceiling(rhx/dcx)-(ncx_o-1-set%nspanlx)
    set%ncovery = ceiling(rhy/dcy)-(ncy_o-1-set%nspanly)
    set%ncoverz = ceiling(rhz/dcz)-(ncz_o-1-set%nspanlz)
    ! Calculate the corners of the CS for checks
    lhcx = dcx*(ncx_o-1-set%nspanlx)
    lhcy = dcy*(ncy_o-1-set%nspanly)
    lhcz = dcz*(ncz_o-1-set%nspanlz)
    rhcx = dcx*(ncx_o-1-set%nspanlx+set%ncoverx)
    rhcy = dcy*(ncy_o-1-set%nspanly+set%ncovery)
    rhcz = dcz*(ncz_o-1-set%nspanlz+set%ncoverz)
    ! First check that the LH and RH corners are large enough
    ! LH
    if(lhx+eps<lhcx) then
       write(*,*) 'CS too small; adjusting nspanlx',lhx,lhcx
       set%nspanlx = set%nspanlx+ceiling((lhcx-lhx)/dcx)
    endif
    if(lhy+eps<lhcy) then
       write(*,*) 'CS too small; adjusting nspanly',lhy,lhcy
       set%nspanly = set%nspanly+ceiling((lhcy-lhy)/dcy)
    endif
    if(lhz+eps<lhcz) then
       write(*,*) 'CS too small; adjusting nspanlz',lhz,lhcz
       set%nspanlz = set%nspanlz+ceiling((lhcz-lhz)/dcz)
    endif
    ! RH
    if(rhx>rhcx+eps) then
       write(*,*) 'CS too small; adjusting ncoverx',rhx,rhcx
       set%ncoverx = set%ncoverx+ceiling((rhx-rhcx)/dcx)
    endif
    if(rhy>rhcy+eps) then
       write(*,*) 'CS too small; adjusting ncovery',rhy,rhcy
       set%ncovery = set%ncovery+ceiling((rhy-rhcy)/dcy)
    endif
    if(rhz>rhcz+eps) then
       write(*,*) 'CS too small; adjusting ncoverz',rhz,rhcz
       set%ncoverz = set%ncoverz+ceiling((rhz-rhcz)/dcz)
    endif
    ! Now check that the CS isn't too large
    ! LH
    if(lhx-lhcx>dcx) then
       write(*,*) 'CS too big; adjusting nspanlx',lhx,lhcx
       set%nspanlx = set%nspanlx-floor((lhx-lhcx)/dcx)
    endif
    if(lhy-lhcy>dcy) then
       write(*,*) 'CS too big; adjusting nspanly',lhy,lhcy
       set%nspanly = set%nspanly-floor((lhy-lhcy)/dcy)
    endif
    if(lhz-lhcz>dcz) then
       write(*,*) 'CS too big; adjusting nspanlz',lhz,lhcz
       set%nspanlz = set%nspanlz-floor((lhz-lhcz)/dcz)
    endif
    ! RH
    if(rhcx-rhx>dcx) then
       write(*,*) 'CS too big; adjusting ncoverx',rhx,rhcx
       set%ncoverx = set%ncoverx-floor((rhx-rhcx)/dcx)
    endif
    if(rhcy-rhy>dcy) then
       write(*,*) 'CS too big; adjusting ncovery',rhy,rhcy
       set%ncovery = set%ncovery-floor((rhy-rhcy)/dcy)
    endif
    if(rhcz-rhz>dcz) then
       write(*,*) 'CS too big; adjusting ncoverz',rhz,rhcz
       set%ncoverz = set%ncoverz-floor((rhz-rhcz)/dcz)
    endif
    return
  end subroutine convert_primary
!!***

!!****f* cover_module/allocate_cs *
!!
!!  NAME 
!!   allocate_cs
!!  USAGE
!! 
!!  PURPOSE
!!   allocates memory to the CS derived type
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   24/11/99
!!  MODIFICATION HISTORY
!!   21/06/2001 dave 
!!    Added ROBODoc header and used cq_abort
!!  SOURCE
!!
  subroutine allocate_cs(set,mx_node,members)

    ! Module usage
    use basic_types
    use GenComms, ONLY: cq_abort
    use global_module, ONLY: area_index
    use memory_module, ONLY: reg_alloc_mem, type_int

    implicit none

    ! Passed variables
    integer :: mx_node
    type(cover_set) :: set
    logical :: members

    ! Local variables
    integer :: stat,irc,ierr

    allocate(set%lab_cell(set%mx_gcover),STAT=stat)
    if(stat/=0) then
       call cq_abort('allocate_cs: error(1)')
    endif
    allocate(set%lab_cover(set%mx_gcover),STAT=stat)
    if(stat/=0) then
       call cq_abort('allocate_cs: error(2)')
    endif
    !TM temporary ???
    allocate(set%inv_lab_cover(set%mx_gcover),STAT=stat)
    if(stat/=0) then
       call cq_abort('allocate_cs: error(3)')
    endif
    !TM temporary ???
    allocate(set%ncover_rem(3*mx_node),STAT=stat)
    if(stat/=0) then
       call cq_abort('allocate_cs: error(4)')
    endif
    call reg_alloc_mem(area_index,3*set%mx_gcover+3*mx_node,type_int)
    if(members) then
       allocate(set%n_ing_cover(set%mx_gcover),STAT=stat)
       if(stat/=0) then
          call cq_abort('allocate_cs: error(5)')
       endif
       allocate(set%icover_ibeg(set%mx_gcover),STAT=stat)
       if(stat/=0) then
          call cq_abort('allocate_cs: error(6)')
       endif
       call reg_alloc_mem(area_index,2*set%mx_gcover,type_int)
       ! These will be allocated later !
       !allocate(set%xcover(set%mx_mcover),STAT=stat)
       !if(stat/=0) then
       !   call cq_abort('allocate_cs: error(7)')
       !endif
       !allocate(set%ycover(set%mx_mcover),STAT=stat)
       !if(stat/=0) then
       !   call cq_abort('allocate_cs: error(8)')
       !endif
       !allocate(set%zcover(set%mx_mcover),STAT=stat)
       !if(stat/=0) then
       !   call cq_abort('allocate_cs: error(9)')
       !endif
    else
       nullify(set%n_ing_cover,set%icover_ibeg, &
            set%xcover,set%ycover,set%zcover)
    endif
    return
  end subroutine allocate_cs
!!***

!!****f* cover_module/deallocate_cs *
!!
!!  NAME 
!!   deallocate_cs
!!  USAGE
!! 
!!  PURPOSE
!!   Deallocates memory from CS derived type
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   24/11/99
!!  MODIFICATION HISTORY
!!   21/06/2001 dave
!!    Added ROBODoc header and cq_abort
!!  SOURCE
!!
  subroutine deallocate_cs(set)

    ! Module usage
    use basic_types
    use GenComms, ONLY: cq_abort
    use global_module, ONLY: area_index
    use memory_module, ONLY: reg_dealloc_mem, type_int

    implicit none

    ! Passed variables
    type(cover_set) :: set

    ! Local variables
    integer :: stat,irc,ierr

    call reg_dealloc_mem(area_index,2*set%mx_gcover,type_int)
    call reg_dealloc_mem(area_index,3*set%mx_gcover,type_int)
    call reg_dealloc_mem(area_index,size(set%ncover_rem),type_int)
    if(ASSOCIATED(set%iprim_group)) then
       call reg_dealloc_mem(area_index,size(set%iprim_group),type_int)
       deallocate(set%iprim_group,STAT=stat)
    endif
    if(stat/=0) then
       call cq_abort('deallocate_cs: error(1)')
    endif
    deallocate(set%ncover_rem,set%zcover,set%ycover,set%xcover, &
         set%lab_cover,set%lab_cell,set%icover_ibeg,set%n_ing_cover, &
         STAT=stat)
    if(stat/=0) then
       call cq_abort('deallocate_cs: error(2)')
    endif
    return
  end subroutine deallocate_cs
!!***

!!****f* cover_module/indexx *
!!
!!  NAME 
!!   indexx - sort
!!  USAGE
!! 
!!  PURPOSE
!!   sbrt indexx: uses the heapsort method to construct a sort-index
!!   array indx() for a supplied integer array iarr(): this means
!!   that the numbers iarr(indx(j)) are in monotonic ascending
!!   order for j=1,2,...nsort. The integer arrays iarr() and indx()
!!   must be dimensioned in the calling routine, and mx_nsort is
!!   the dimension. A test is made to ensure that n.le.mx_nsort, and if
!!   this fails, sbrt indexx stops with an error message. This routine
!!   is taken, with slight modifications, from Numerical Recipes,
!!   1st edition, p. 233.
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   24/11/99
!!  MODIFICATION HISTORY
!!   21/06/2001 dave
!!    Added ROBODoc header and used cq_abort
!!  SOURCE
!!
  subroutine indexx(mx_nsort,nsort,iarr,indx)

    ! Module usage
    use datatypes
    use GenComms, ONLY: cq_abort

    implicit none

    ! Passed variables
    integer mx_nsort, nsort
    integer :: iarr(mx_nsort),indx(mx_nsort)

    ! Local variables
    integer :: i,j,m,ir,indxt,iq,irc,ierr

    ! check on array bounds 
    if(nsort>mx_nsort.OR.nsort<1) then
       call cq_abort('indexx: too many elements to be sorted',nsort,mx_nsort)
    endif
    do j=1,nsort
       indx(j)=j
    enddo
    if(nsort.eq.1) return
    m=1+nsort/2
    ir=nsort
    do while(.true.)
       if(m>1) then
          m=m-1
          indxt=indx(m)
          iq=iarr(indxt)
       else
          if(ir<1.OR.ir>mx_nsort) call cq_abort("Overflow in indexx: ",ir,mx_nsort)
          indxt=indx(ir)
          iq=iarr(indxt)
          indx(ir)=indx(1)
          ir=ir-1
          if(ir==1) then
             indx(1)=indxt
             return
          endif
       endif
       i=m
       j=m+m
       do while(j<=ir) 
          if(j<ir) then
             if(iarr(indx(j))<iarr(indx(j+1))) j=j+1
          endif
          if(iq<iarr(indx(j))) then
             indx(i)=indx(j)
             i=j
             j=j+j
          else
             j=ir+1
          endif
       enddo
       indx(i)=indxt
    enddo
  end subroutine indexx
!!***
end module cover_module
