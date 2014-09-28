! $Id: UpdateMember_module.f90 81 2013-07-05 michi $
! -----------------------------------------------------------
! Module UpdateMember_module
! -----------------------------------------------------------
! Code area 8: indexing
! -----------------------------------------------------------

!!****h* Conquest/group_module
!!  NAME
!!   UpdateMember_module
!!  PURPOSE
!!   Updates all the necessary data (variables) when running MD
!!   or structure optimisation
!!  AUTHOR
!!   Michiaki Arita
!!  CREATION DATE
!!   2013/07/02
!!  MODIFICATION HISTORY
!!   2014/09/15 18:30 lat
!!    fixed call start/stop_timer to timer_module (not timer_stdlocks_module !)
!!  SOURCE
!!
module UpdateMember_module

  ! Module usage
  use datatypes
  use basic_types
  use global_module,          only: flag_MDdebug, iprint_MDdebug
  use timer_module,           only: start_timer, stop_timer
  use timer_stdclocks_module, only: tmr_std_indexing, tmr_std_allocation

  implicit none
  save

  type(group_set) :: parts  ! Partitions of atoms
  type(group_set) :: blocks ! Blocks of integration grid points

  character(80),private :: RCSid = "$Id$"

!!***

contains

  !!****f* UpdateMember_module/group_update_mparts *
  !!
  !!  NAME 
  !!   group_update_mparts
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Updates the necessary information on group_set 'parts'.
  !!   All the updates of primary_set and cover_set are besed on 
  !!   this process.
  !!  INPUTS
  !!   velocity
  !!     - When running MD, velocity is precisely particle velocities, 
  !!       but in the case of 'cg', velocity is direction(3,ni_in_cell)
  !!       ,which is a searching direction for conjugate gradient.
  !!   iteration (optional)
  !!     - We have iteration only when running MD. The code writes out 
  !!       the current MD step in 'make_prt.dat'.
  !!     - [NOTE:] This may be deleted in the next update because showing
  !!               the MD step in 'make_prt.dat' is NOT a good idea - 
  !!               it should be written out, e.g. InfoGlobal.dat .
  !!  USES
  !!   basic_types,global_module,species_module,atom_dispenser,group_module
  !!  AUTHOR
  !!   Michiaki Arita
  !!  CREATION DATE
  !!   2013/07/02
  !!  MODIFICATION HISTORY
  !!   2013/08/21 M.Arita
  !!   - Removed iteration
  !!  SOURCE
  !!
  !ORI subroutine group_update_mparts(velocity,iteration)
  subroutine group_update_mparts(velocity)

    ! Module usage
    use basic_types
    use GenComms, ONLY: cq_abort
    use global_module, ONLY: ni_in_cell,id_glob,id_glob_inv, &
                             id_glob_old,id_glob_inv_old, species_glob, &
                             x_atom_cell,y_atom_cell,z_atom_cell,atom_coord
    use species_module, ONLY: species
    use atom_dispenser, ONLY: allatom2part
    use group_module, ONLY: parts
    ! DB
    use io_module, ONLY: get_file_name
    use input_module, ONLY: io_assign,io_close
    use GenComms, ONLY: ionode,inode
    use global_module, ONLY: numprocs
    use maxima_module, ONLY: maxpartsproc,maxatomspart

    implicit none

    ! passed variables
    !ORI integer,intent(in), optional :: iteration
    real(double) :: velocity(3,ni_in_cell)

    ! local variables
    integer :: mx_mem_grp,np_cc,np_cc_old,CC_part,np_old,ia_part_lab
    integer :: np,ni,stat_alloc,iglob,ni_old,natom,ind_group,nnd
    integer, allocatable :: ind_part(:),id_atom_old(:)
    integer, allocatable :: ind_ia(:),ind_ip(:)
    real(double), allocatable :: velocity_tmp(:,:)
    real(double), allocatable :: x_atom_cell_tmp(:)
    real(double), allocatable :: y_atom_cell_tmp(:)
    real(double), allocatable :: z_atom_cell_tmp(:)
    logical :: flag_init,flag
    logical, allocatable :: flag_hitprts(:)
    character(12) :: file_prt = 'make_prt.dat'
    integer :: lun_prt,inode_beg

    ! DB
    integer :: lun_db,stat,id_global
    character(20) :: file_name

    !! ---------- DEBUG ---------- !!
    ! NOTE: This file will be sizable.
    if (flag_MDdebug) then
      call get_file_name('gp_update',numprocs,inode,file_name)
      call io_assign(lun_db)
      open  (lun_db,file=file_name,position='append')
      write (lun_db,*) "--- Before Update ---"
      write (lun_db,*) "parts%mx_mem_grp:", parts%mx_mem_grp
      write (lun_db,*) "parts%nm_group :", parts%nm_group(1:)
      write (lun_db,*) "parts%icell_beg:", parts%icell_beg(1:)
      write (lun_db,*) "id_glob        :", id_glob(1:)
      write (lun_db,*) "id_glob_old    :", id_glob_old(1:)
      write (lun_db,*) "id_glob_inv    :", id_glob_inv(1:)
      write (lun_db,*) "id_glob_inv_old:", id_glob_inv_old(1:)
      write (lun_db,*) ""
      write (lun_db,*) "velocity -- Before updating."
      do ni = 1, ni_in_cell
        write (lun_db,*) ni, velocity(1,ni), velocity(2,ni), velocity(3,ni)
      enddo
      write (lun_db,*) "Atomic positions -- Before updating."
      do ni = 1, ni_in_cell
        write (lun_db,*) ni, x_atom_cell(ni), y_atom_cell(ni), z_atom_cell(ni)
      enddo
    endif
    !! ---------- DEBUG ---------- !!

    ! NOTE: 
    ! velocity(1:3, ni_in_cell) stand for exactly velocities in MD and "direction(1:3,ni_in_cell)" in CG !
    allocate (id_atom_old(ni_in_cell),ind_ia(ni_in_cell),ind_ip(ni_in_cell), STAT=stat_alloc)
    if (stat_alloc.NE.0) call cq_abort('Error allocating id_atom_old:', ni_in_cell)
    allocate (flag_hitprts(parts%mx_gcell), STAT=stat_alloc)
    if (stat_alloc.NE.0) call cq_abort('Error: allocating flag_hitprts:', parts%mx_gcell)

    ! Hold the old relationship between global and partition labellings.
    id_glob_old = id_glob
    id_glob_inv_old = id_glob_inv

    ! Get the partitions in sim-cell (CC) for all atoms (glob).
    allocate (ind_part(ni_in_cell), STAT=stat_alloc)
    if (stat_alloc.NE.0) call cq_abort('Error allocating ind_part: ', ni_in_cell)
    call allatom2part(ind_part) ! --> parts. in sim cell in CC.

    parts%nm_group = 0
    ind_ia = 0
    ind_ip = 0
    do iglob = 1, ni_in_cell                                            ! iglob: global labelling.
      CC_part = ind_part(iglob)
      parts%nm_group(CC_part) = parts%nm_group(CC_part) + 1
      ind_ia(iglob) = parts%nm_group(CC_part)                           ! atomID in partition "CC_part" 
                                                                        !   -> used to get updated id_glob(:)
      if (ind_ia(iglob).LE.0 .OR. ind_ia(iglob).GT.ni_in_cell) &
        call cq_abort('Error in ind_ia(iglob)', iglob, ind_ia(iglob))
      ind_ip(iglob) = CC_part                                           ! parts in sim-cell (CC).
      if (ind_ip(iglob).LE.0 .OR. ind_ip(iglob).GT.parts%mx_gcell) &
        call cq_abort('Error in ind_ip(iglob)', iglob, ind_ip(iglob))
    enddo
    parts%mx_mem_grp = maxval(parts%nm_group)
    maxatomspart = parts%mx_mem_grp

    !! ------------ DEBG: ------------ !!
    if (flag_MDdebug) then
      write (lun_db, '(/a)') "ind_ia(1:ni_in_cell):"
      write (lun_db,*) ind_ia(1:)
      write (lun_db, '(a)' ) "ind_ip(1:ni_in_cell):"
      write (lun_db,*) ind_ip(1:)
      write (lun_db,*) ""
    endif
    !! ------------ DEBG: ------------ !!

    ! Get icell_beg(:).
    parts%icell_beg = 0
    flag_init = .true.
    parts%icell_beg=1 ; np_old=1
    flag=.false.
    do nnd = 1, numprocs
      if (parts%ng_on_node(nnd).GT.0) then
        do np = 1, parts%ng_on_node(nnd)
          CC_part=parts%ngnode(parts%inode_beg(nnd)+np-1)
          if (.NOT. flag) then
            parts%icell_beg(CC_part)=parts%icell_beg(np_old)
          else
            parts%icell_beg(CC_part)=parts%icell_beg(np_old)+parts%nm_group(np_old)
          endif
          ! No atoms in part?
          if (parts%nm_group(CC_part).GT.0) then
            flag=.true.
          else
            flag=.false.
          endif
          np_old=CC_part
          ! Error check.
          !if (parts%icell_beg(CC_part).GT.ni_in_cell) &
          !call cq_abort('Error: icell_beg must not be larger than ni_in_cell:', &
          !               parts%icell_beg(CC_part))
        enddo
      endif
    enddo

    ! Update 'id_glob', 'id_glob_inv', 'species', 'xyz_atom_cell' and velocity/direction.
    do iglob = 1, ni_in_cell
      ia_part_lab = parts%icell_beg(ind_ip(iglob))+ind_ia(iglob)-1 ! parts. labelling
      if (ia_part_lab .LE.0 .OR. ia_part_lab .GT. ni_in_cell) &
        call cq_abort('Error in ia_part_lab:', ia_part_lab)
      id_glob(ia_part_lab) = iglob
      id_glob_inv(iglob)   = ia_part_lab
    enddo
    allocate (velocity_tmp(3,ni_in_cell), STAT=stat_alloc)
    if (stat_alloc.NE.0) call cq_abort('Error allocating velocity_tmp:', 3, ni_in_cell)
    allocate (x_atom_cell_tmp(ni_in_cell),y_atom_cell_tmp(ni_in_cell),z_atom_cell_tmp(ni_in_cell), &
              STAT=stat_alloc)
    if (stat_alloc.NE.0) call cq_abort('Error allocating xyz_atom_cell_tmp:', ni_in_cell)
    ! Hold old relations
    velocity_tmp    = velocity
    x_atom_cell_tmp = x_atom_cell
    y_atom_cell_tmp = y_atom_cell
    z_atom_cell_tmp = z_atom_cell
    do ni = 1, ni_in_cell                                          ! parts. labelling
      id_global = id_glob(ni)
      ni_old = id_glob_inv_old(id_global)
      if (flag_MDdebug) then
        write (lun_db,*) "globID et ni_old:", id_global, ni_old                   !DB
        write (lun_db,*) velocity(1,ni_old),velocity(2,ni_old),velocity(3,ni_old) !DB
      endif
      x_atom_cell(ni) = x_atom_cell_tmp(ni_old)                    ! Unwrapped
      y_atom_cell(ni) = y_atom_cell_tmp(ni_old)                    ! Unwrapped
      z_atom_cell(ni) = z_atom_cell_tmp(ni_old)                    ! Unwrapped
      ! The following updates are already done in velocityVerlet -- See update_atom_coord.
      !NOT NEEDED atom_coord(1,id_global) = x_atom_cell(ni)
      !NOT NEEDED atom_coord(2,id_global) = y_atom_cell(ni)
      !NOT NEEDED atom_coord(3,id_global) = z_atom_cell(ni)
      velocity(1,ni)  = velocity_tmp(1,ni_old)
      velocity(2,ni)  = velocity_tmp(2,ni_old)
      velocity(3,ni)  = velocity_tmp(3,ni_old)
      species(ni)     = species_glob(id_global)
    enddo

    ! Check whether bundle is not empty.
    ! This might be critical in some cases.
    !  --> Needs to redistribute atoms to processors.
    if (inode.EQ.ionode) then
      do nnd = 1, numprocs
        natom = 0
        if (parts%ng_on_node(nnd).GT.0) then
          do np = 1, parts%ng_on_node(nnd)
            ind_group = parts%ngnode(parts%inode_beg(nnd)+np-1)  ! CC
            natom = natom + parts%nm_group(ind_group)
          enddo
          if (natom.EQ.0) &
            call cq_abort('Critical error: No atoms in bundle!', nnd)
        endif
      enddo !(nnd, numprocs)
    endif

    deallocate (ind_part, STAT=stat_alloc)
    if (stat_alloc.NE.0) call cq_abort('Error deallocating ind_part: ', ni_in_cell)
    deallocate (id_atom_old,ind_ia,ind_ip, STAT=stat_alloc)
    if (stat_alloc.NE.0) call cq_abort('Error deallocating id_atom_old: ', ni_in_cell)
    deallocate (flag_hitprts, STAT=stat_alloc)
    if (stat_alloc.NE.0) call cq_abort('Error deallocating flag_hitprts: ', parts%mx_gcell)
    deallocate (velocity_tmp, STAT=stat_alloc)
    if (stat_alloc.NE.0) call cq_abort('Error deallocating velocity_tmp:', 3, ni_in_cell)
    deallocate (x_atom_cell_tmp,y_atom_cell_tmp,z_atom_cell_tmp, STAT=stat_alloc)
    if (stat_alloc.NE.0) call cq_abort('Error deallocating xyz_atom_cell_tmp:', ni_in_cell)

    !! ---------- DEBUG ---------- !!
    if (flag_MDdebug) then
      write (lun_db,*) ""
      write (lun_db,*) "--- After Update ---"
      write (lun_db,*) "parts%mx_mem_grp:", parts%mx_mem_grp
      write (lun_db,*) "parts%nm_group  :", parts%nm_group(1:)
      write (lun_db,*) "parts%icell_beg :", parts%icell_beg(1:)
      write (lun_db,*) "id_glob         :", id_glob(1:)
      write (lun_db,*) "id_glob_old     :", id_glob_old(1:)
      write (lun_db,*) "id_glob_inv     :", id_glob_inv(1:)
      write (lun_db,*) "id_glob_inv_old :", id_glob_inv_old(1:)
      do np = 1, parts%mx_gcell
        if (parts%nm_group(np).GT.0) then
          do ni = 1, parts%nm_group(np)
            id_global = id_glob(parts%icell_beg(np)+ni-1)
            write (lun_db,*) "CCprts, id_glob, ni:", np,id_global, &
                                                     parts%icell_beg(np)+ni-1
          enddo
        endif
      enddo
      write (lun_db,*) ""
      write (lun_db,*) "velocity -- After updating."
      write (lun_db,*) "ni, id_glob, velocity(1:3, ni)"
      do ni = 1, ni_in_cell
        !write (lun_db,*) ni, id_glob(ni), velocity(1,ni), velocity(2,ni), velocity(3,ni)
        write (lun_db,*) velocity(1,ni), velocity(2,ni), velocity(3,ni)
      enddo
      write (lun_db,*) "Atomic positions -- After updating."
      do ni = 1, ni_in_cell
        write (lun_db,*) ni, x_atom_cell(ni), y_atom_cell(ni), z_atom_cell(ni)
      enddo
      write (lun_db,*) ""
      call io_close(lun_db)
    endif
    !! ---------- DEBUG ---------- !!

    ! Output 'make_prt.dat'
    if (inode.EQ.ionode) then
      call io_assign(lun_prt)
      open (lun_prt,file=file_prt,iostat=stat)
      if (stat.NE.0) call cq_abort('Error opening make_prt.dat in group_update_mparts!')
      write (lun_prt,1) parts%ngcellx,parts%ngcelly,parts%ngcellz
      write (lun_prt,2) numprocs
      inode_beg=1
      do nnd = 1, numprocs
        write (lun_prt,3) nnd,parts%ng_on_node(nnd),inode_beg
        inode_beg=inode_beg+parts%ng_on_node(nnd)
        if (parts%ng_on_node(nnd).GT.0) then
          do np = 1, parts%ng_on_node(nnd)
            ! The following used to be wrong. Now ok.
            CC_part=parts%ngnode(parts%inode_beg(nnd)+np-1)
            write (lun_prt,4) np,CC_part,parts%nm_group(CC_part),parts%icell_beg(CC_part)
            if (parts%nm_group(CC_part).GT.0) then
              do ni = 1, parts%nm_group(CC_part)
                id_global=id_glob(parts%icell_beg(CC_part)+ni-1)
                write (lun_prt,5) ni,id_global
              enddo
            endif
          enddo
        endif
      enddo
      ! This may be deleted in the next update
      !ORI if (present(iteration)) write (lun_prt,*) iteration
      call io_close(lun_prt)
    endif

1   format(3i5)
2   format(i10)
3   format(3i10)
4   format(5x,4i10)
5   format(2i8)

    return
  end subroutine group_update_mparts
  !!***

  !!****f* UpdateMember_module/deallocate_PSmember *
  !!
  !!  NAME 
  !!   deallocate_PSmember
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Deallocates the necessary arrays of 'bundle' for updates
  !!  INPUTS
  !! 
  !!  USES
  !!   primary_module
  !!  AUTHOR
  !!   Michiaki Arita
  !!  CREATION DATE
  !!   2013/07/02
  !!  MODIFICATION HISTORY
  !!
  !!  SOURCE
  !!
  subroutine deallocate_PSmember

    ! Module usage
    use basic_types
    use GenComms, ONLY: cq_abort
    use primary_module, ONLY: bundle
    ! DB
    use io_module, ONLY: get_file_name
    use GenComms, ONLY: inode
    use global_module, ONLY: numprocs
    use input_module, ONLY: io_assign,io_close

    implicit none

    ! local variables
    integer :: stat_alloc

    ! DB
    integer :: lun_db
    character(20) :: file_name

    !! ---------- DEBUG ---------- !!
    if (flag_MDdebug .AND. iprint_MDdebug.GT.3) then
      call get_file_name('deallocPS',numprocs,inode,file_name)
      call io_assign(lun_db)
      open (lun_db,file=file_name)
      write (lun_db,*) "--- Befor update. ---"
      write (lun_db,*) "mx_iprim, n_prim:", bundle%mx_iprim, bundle%n_prim
      write (lun_db,*) "mx_ngonn   :", bundle%mx_ngonn
      write (lun_db,*) "nm_nodgroup:", bundle%nm_nodgroup(1:)
      write (lun_db,*) "nm_nodbeg  :", bundle%nm_nodbeg(1:)
      write (lun_db,*) "iprim_seq  :", bundle%iprim_seq(1:)
      write (lun_db,*) "ig_prim    :", bundle%ig_prim(1:)
      write (lun_db,*) "xprim      :", bundle%xprim(1:)
      write (lun_db,*) "yprim      :", bundle%yprim(1:)
      write (lun_db,*) "zprim      :", bundle%zprim(1:)
      write (lun_db,*) "species    :", bundle%species(1:)
      call io_close(lun_db)
    endif
    !! ---------- DEBUG ---------- !!

    deallocate(bundle%ig_prim,bundle%xprim,bundle%yprim,bundle%zprim, &
               bundle%iprim_seq,bundle%species, STAT=stat_alloc) 
    if(stat_alloc.NE.0) &
       call cq_abort('Error deallocating bundle associated with mebers.')

    return
  end subroutine deallocate_PSmember
  !!***

  !!****f* UpdateMember_module/allocate_PSmember *
  !!
  !!  NAME 
  !!   allocate_PSmember
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Allocates all the necessary arrays for bundle
  !!  INPUTS
  !!   groups
  !!     - 'parts' is always called since this subroutine updates
  !!       all the necessary information related to 'bundle'
  !!  USES
  !!   primary_module,maxima_module
  !!  AUTHOR
  !!   Michiaki Arita
  !!  CREATION DATE
  !!   2013/07/02
  !!  MODIFICATION HISTORY
  !!
  !!  SOURCE
  !!
  subroutine allocate_PSmember(groups)

    ! Module usage
    use basic_types
    use global_module, ONLY: ni_in_cell
    use GenComms, ONLY: cq_abort,inode,gmax,my_barrier
    use primary_module, ONLY: bundle
    use maxima_module, ONLY: maxatomsproc
    ! DB
    use io_module, ONLY: get_file_name
    use global_module, ONLY: numprocs
    use input_module, ONLY: io_assign,io_close

    implicit none

    ! passed variables
    type(group_set), target :: groups

    ! local variables
    integer :: ng,stat_alloc
    integer :: ind_group,n_prim,mx_iprim

    ! DB
    integer :: lun_db
    character(20) :: file_name

    bundle%groups => groups

    !! ---------- DEBUG ---------- !!
    if (flag_MDdebug .AND. iprint_MDdebug.GT.3) then
      call get_file_name('allocate_PS',numprocs,inode,file_name)
      call io_assign(lun_db)
      open (lun_db, file=file_name)
      write (lun_db,*) "Check if parts is updated and read properly in bundle!"
      write (lun_db,*) "parts%mx_mem_grp :", groups%mx_mem_grp
      write (lun_db,*) "parts%nm_group   :", groups%nm_group(1:)
      write (lun_db,*) "parts%icell_beg  :", groups%icell_beg(1:)
      write (lun_db,*) ""
      write (lun_db,*) ""
      write (lun_db,*) "--- Before update."
      write (lun_db,*) "mx_iprim    :", bundle%mx_iprim
      write (lun_db,*) "n_prim      :", bundle%n_prim
      write (lun_db,*) "maxatomsproc:", maxatomsproc
    endif
    !! ---------- DEBUG ---------- !!

    n_prim = 0
    do ng = 1, groups%ng_on_node(inode)
      ind_group = groups%ngnode(groups%inode_beg(inode)+ng-1)  ! CC in sim-cell.
      n_prim = n_prim + groups%nm_group(ind_group)
    enddo
    bundle%n_prim = n_prim
    call my_barrier()
    call gmax(n_prim)
    mx_iprim = n_prim
    bundle%mx_iprim = mx_iprim

    allocate (bundle%iprim_seq(mx_iprim),bundle%ig_prim(mx_iprim), &
              bundle%xprim(mx_iprim),bundle%yprim(mx_iprim),       &
              bundle%zprim(mx_iprim),bundle%species(mx_iprim),     &
              STAT=stat_alloc)
    if (stat_alloc.NE.0) &
      call cq_abort('Error allocating bundle associated with members: ', mx_iprim)

    ! This should be bundle%n_prim...
    ! NOTE: In the old version of CONQUEST, maxatomsproc was constant
    !       during a job.
    maxatomsproc = bundle%mx_iprim

    !! ---------- DEBUG --------- !!
    if (flag_MDdebug .AND. iprint_MDdebug.GT.3) then
      write (lun_db,*) ""
      write (lun_db,*) "--- After update."
      write (lun_db,*) "mx_iprim        : ", bundle%mx_iprim
      write (lun_db,*) "n_prim          : ", bundle%n_prim
      write (lun_db,*) "maxatomsproc    : ", maxatomsproc
      write (lun_db,*) "bundle%nx_origin: ", bundle%nx_origin
      write (lun_db,*) "bundle%ny_origin: ", bundle%ny_origin
      write (lun_db,*) "bundle%nz_origin: ", bundle%nz_origin
      call io_close(lun_db)
    endif
    !! ---------- DEBUG --------- !!

    return
  end subroutine allocate_PSmember
  !!***

  !!****f* UpdateMember_module/primary_update_mparts *
  !!
  !!  NAME 
  !!   primary_update_mparts
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Updates the member information in bundle
  !!  INPUTS
  !! 
  !!  USES
  !!   global_module,group_module,primary_module,species_module
  !!  AUTHOR
  !!   Michiaki Arita
  !!  CREATION DATE
  !!   2013/07/02
  !!  MODIFICATION HISTORY
  !!
  !!  SOURCE
  !!
  subroutine primary_update_mparts

    ! Module usage
    use basic_types
    use global_module, ONLY: rcellx,rcelly,rcellz,x_atom_cell,y_atom_cell,z_atom_cell, &
                             id_glob
    use GenComms, ONLY: cq_abort,inode
    use group_module, ONLY: parts
    use primary_module, ONLY: bundle
    use species_module, ONLY: species
    ! DB
    use io_module, ONLY: get_file_name
    use global_module, ONLY: numprocs,ni_in_cell
    use input_module, ONLY: io_assign,io_close

    implicit none

    ! local variables
    integer :: ng,ni,stat_alloc
    integer :: ind_group,nx,ny,nz,nx1,ny1,nz1,n_prim
    real(double) :: dcellx,dcelly,dcellz,xadd,yadd,zadd

    ! DB
    integer :: lun_db
    character(20) :: file_name

    !! ---------- DEBUG ---------- !!
    if (flag_MDdebug .AND. iprint_MDdebug.GT.3) then
      call get_file_name('updatePS',numprocs,inode,file_name)
      call io_assign(lun_db)
      open (lun_db,file=file_name)
    endif
    !! ---------- DEBUG ---------- !!

    dcellx = rcellx / parts%ngcellx
    dcelly = rcelly / parts%ngcelly
    dcellz = rcellz / parts%ngcellz
    bundle%nm_nodgroup=0  ; bundle%nm_nodbeg=0
    bundle%nm_nodbeg(1)=1 ; n_prim = 0
    do ng = 1, parts%ng_on_node(inode)
      ind_group = parts%ngnode(parts%inode_beg(inode)+ng-1) ! CC in sim-cell.
      ! Decompose 'ind_group' into offset.
      nx=1+(ind_group-1)/(parts%ngcelly*parts%ngcellz)
      ny=1+(ind_group-1-(nx-1)*parts%ngcelly*parts%ngcellz)/parts%ngcellz
      nz=ind_group-(nx-1)*parts%ngcelly*parts%ngcellz-(ny-1)*parts%ngcellz
      ! Offset of 'ng' w.r.t. Origin of PS.
      nx1=bundle%nx_origin+bundle%idisp_primx(ng)
      ny1=bundle%ny_origin+bundle%idisp_primy(ng)
      nz1=bundle%nz_origin+bundle%idisp_primz(ng)
      ! Offset to bring atoms to PS.
      xadd = real(nx1-nx,double)*dcellx
      yadd = real(ny1-ny,double)*dcelly
      zadd = real(nz1-nz,double)*dcellz
      bundle%nm_nodgroup(ng) = parts%nm_group(ind_group)
      if (bundle%nm_nodgroup(ng).GT.0) then
        do ni = 1, bundle%nm_nodgroup(ng)
          n_prim = n_prim + 1
          if (n_prim.GT.bundle%mx_iprim) &
            call cq_abort('Error: n_prim exceeds mx_iprim.')
          bundle%iprim_seq(n_prim) = ni
          bundle%ig_prim(n_prim) = id_glob(parts%icell_beg(ind_group)+ni-1)
          bundle%xprim(n_prim)   = x_atom_cell(parts%icell_beg(ind_group)+ni-1) + xadd
          bundle%yprim(n_prim)   = y_atom_cell(parts%icell_beg(ind_group)+ni-1) + yadd
          bundle%zprim(n_prim)   = z_atom_cell(parts%icell_beg(ind_group)+ni-1) + zadd
          bundle%species(n_prim) = species(parts%icell_beg(ind_group)+ni-1)
          if (flag_MDdebug .AND. iprint_MDdebug.GT.3) &
            write (lun_db,*) "icell_beg+ni-1:", parts%icell_beg(ind_group)+ni-1 !DB
        enddo
      endif
      if (ng.LT.parts%ng_on_node(inode)) &
        bundle%nm_nodbeg(ng+1) = bundle%nm_nodbeg(ng) + bundle%nm_nodgroup(ng)
    enddo

    !! ---------- DEBUG ---------- !!
    if (flag_MDdebug .AND. iprint_MDdebug.GT.3) then
      write (lun_db,*) ""
      write (lun_db,*) "xyz_atom_cell:"
      do ni = 1, ni_in_cell
        write (lun_db,*) x_atom_cell(ni),y_atom_cell(ni),z_atom_cell(ni)
      enddo
      write (lun_db,*) "--- After update. ---"
      write (lun_db,*) "mx_iprim, n_prim:", bundle%mx_iprim, bundle%n_prim
      write (lun_db,*) "mx_ngonn   :", bundle%mx_ngonn
      write (lun_db,*) "nm_nodgroup:", bundle%nm_nodgroup(1:)
      write (lun_db,*) "nm_nodbeg  :", bundle%nm_nodbeg(1:)
      write (lun_db,*) "iprim_seq  :", bundle%iprim_seq(1:)
      write (lun_db,*) "ig_prim    :", bundle%ig_prim(1:)
      write (lun_db,*) "xprim      :", bundle%xprim(1:)
      write (lun_db,*) "yprim      :", bundle%yprim(1:)
      write (lun_db,*) "zprim      :", bundle%zprim(1:)
      write (lun_db,*) "species    :", bundle%species(1:)
      !write (lun_db,*) "parts%nm_group   :", bundle%parts%nm_group(1:)
      !write (lun_db,*) "parts%icell_beg  :", bundle%parts%icell_beg(1:)
      call io_close(lun_db)
    endif
    !! ---------- DEBUG ---------- !!

    return
  end subroutine primary_update_mparts
  !!***

  !!****f* UpdateMember_module/deallocate_CSmember *
  !!
  !!  NAME 
  !!   deallocate_CSmember
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Deallocates all the necessary arrays on covering sets for 
  !!   member-updates
  !!  INPUTS
  !!   set
  !!    - BCS_parts,DCS_parts,ewald_CS,D2_CS
  !!  USES
  !!  
  !!  AUTHOR
  !!   Michiaki Arita
  !!  CREATION DATE
  !!   2013/07/02
  !!  MODIFICATION HISTORY
  !!
  !!  SOURCE
  !!
  subroutine deallocate_CSmember(set)
  !ORI subroutine deallocate_CSmember(set,prim)

    ! Module usage
    use basic_types
    use GenComms, ONLY: cq_abort
    ! DB
    use io_module, ONLY: get_file_name
    use GenComms, ONLY: inode
    use global_module, ONLY: numprocs
    use input_module, ONLY: io_assign,io_close

    implicit none

    ! passed variables
    type(cover_set)   :: set

    ! local variables
    integer :: stat_alloc

    ! DB
    integer :: lun_db
    character(20) :: file_name

    deallocate(set%xcover, STAT=stat_alloc)
    if(stat_alloc.NE.0) &
       call cq_abort('Error deallocating set%xcover:', set%mx_mcover)
    deallocate(set%ycover, STAT=stat_alloc)
    if(stat_alloc.NE.0) &
       call cq_abort('Error deallocating set%ycover:', set%mx_mcover)
    deallocate(set%zcover, STAT=stat_alloc)
    if(stat_alloc.NE.0) &
       call cq_abort('Error deallocating set%zcover:', set%mx_mcover)
    deallocate(set%ig_cover, STAT=stat_alloc)
    if(stat_alloc.NE.0) &
       call cq_abort('Error deallocating set%ig_cover:', set%mx_mcover)
    if (associated(set%iprim_group)) then
      deallocate(set%iprim_group, STAT=stat_alloc)
      if(stat_alloc.NE.0) &
         !ORI call cq_abort('Error deallocating set%iprim_group:', prim%mx_iprim)
         call cq_abort('Error deallocating set%iprim_group:')
    endif

    return
  end subroutine deallocate_CSmember
  !!***

  !!****f* UpdateMember_module/allocate_CSmember *
  !!
  !!  NAME 
  !!   allocate_CSmember
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Allocates all the necessary arrays on covering sets
  !!   for member-updates
  !!  INPUTS
  !!   set,groups,nx_in_cover,ny_in_cover,nz_in_cover,nmodx,mody,nmodz,
  !!   dcellx,dcelly,dcellz,prim
  !!  USES
  !!   basic_types,global_module,cover_module
  !!  AUTHOR
  !!   Michiaki Arita
  !!  CREATION DATE
  !!    2013/07/02
  !!  MODIFICATION HISTORY
  !!
  !!  SOURCE
  !!
  subroutine allocate_CSmember(set,groups,nx_in_cover,ny_in_cover,nz_in_cover, &
                               nmodx,nmody,nmodz,dcellx,dcelly,dcellz,prim)

    ! Module usage
    use basic_types
    use global_module, ONLY: rcellx,rcelly,rcellz
    use GenComms, ONLY: cq_abort
    use cover_module, ONLY: indexx
    ! DB
    use io_module, ONLY: get_file_name
    use GenComms, ONLY: inode
    use global_module, ONLY: numprocs
    use input_module, ONLY: io_assign,io_close

    implicit none

    ! passed variables
    integer :: nmodx,nmody,nmodz
    ! x,y,z numbering of CS groups (CC)
    integer, allocatable :: nx_in_cover(:),ny_in_cover(:),nz_in_cover(:)
    type(cover_set)          :: set
    type(primary_set),target,optional :: prim
    type(group_set),target   :: groups

    ! local variables
    ! No. of periodic images of a group.
    integer, allocatable :: nrepx(:),nrepy(:),nrepz(:)

    ! variables for irreducible CS.
    integer, allocatable :: ind_min(:)
    integer, allocatable :: ngcx_min(:)
    integer, allocatable :: ngcy_min(:)
    integer, allocatable :: ngcz_min(:)
    integer, allocatable :: min_sort(:)

    integer :: noccx,noccy,noccz,nremx,nremy,nremz,minx,miny,minz,ngcx,ngcy,ngcz
    integer :: ng_in_cell,ng_in_min,ind,nqx,nqy,nqz,ind_qart,ino
    integer :: ind_cover,nrx,nry,nrz
    integer :: nm_in_cover,nsx,nsy,nsz
    integer :: stat_alloc
    real(double) :: dcellx,dcelly,dcellz

    ! DB
    integer :: lun_db
    character(20) :: file_name

    ! Allocate set%iprim_group.
    set%groups => groups
    set%prim   => prim
    if (present(prim)) then
      allocate(set%iprim_group(prim%mx_iprim), STAT=stat_alloc)
      if (stat_alloc.NE.0) &
        call cq_abort('Error allocating iprim_group: ', prim%mx_iprim)
    endif

    ! Get set%mx_mcover.
    ! Allocation.
    allocate (nx_in_cover(set%ng_cover),ny_in_cover(set%ng_cover), &
              nz_in_cover(set%ng_cover), STAT=stat_alloc)
    if (stat_alloc.NE.0) & 
      call cq_abort('Error allocating nx_in_cover:', set%ng_cover)
    allocate (nrepx(groups%mx_gedge),nrepy(groups%mx_gedge),nrepz(groups%mx_gedge), &
              STAT=stat_alloc)
    if (stat_alloc.NE.0) &
      call cq_abort('Error allocating nrepx,y,z:', groups%mx_gedge)
    allocate (ind_min(groups%mx_gcell),ngcx_min(groups%mx_gcell),ngcy_min(groups%mx_gcell), &
              ngcz_min(groups%mx_gcell),min_sort(groups%mx_gcell), STAT=stat_alloc)
    if (stat_alloc.NE.0) &
      call cq_abort('Error allocating ind_min,ngcx,y,z_min,min_sort:', groups%mx_gcell)

    ! Conversion factors from unit cell lengths -> groups.
    dcellx = rcellx/real(groups%ngcellx,double)
    dcelly = rcelly/real(groups%ngcelly,double)
    dcellz = rcellz/real(groups%ngcellz,double)

    ! Refer to make_cs.
    nmodx = ((groups%ngcellx+set%nspanlx-1)/groups%ngcellx)*groups%ngcellx
    nmody = ((groups%ngcelly+set%nspanly-1)/groups%ngcelly)*groups%ngcelly
    nmodz = ((groups%ngcellz+set%nspanlz-1)/groups%ngcellz)*groups%ngcellz

    ! Create apparatus for p.b.c. .
    ! X-direction.
    noccx = set%ncoverx / groups%ngcellx
    nremx = set%ncoverx - noccx*groups%ngcellx
    minx  = min(set%ncoverx,groups%ngcellx)
    if (minx.GT.groups%mx_gedge) &
      call cq_abort('allocate_CS: too many groups in x-edgs', minx)
    do ngcx = 1, minx
      if (ngcx.LE.nremx) then
        nrepx(ngcx) = noccx + 1
      else
        nrepx(ngcx) = noccx
      endif
    enddo
    ! Y-direction.
    noccy = set%ncovery / groups%ngcelly
    nremy = set%ncovery - noccy*groups%ngcelly
    miny  = min(set%ncovery,groups%ngcelly)
    if (miny.GT.groups%mx_gedge) &
      call cq_abort('allocate_CS: too many groups in y-edgs', miny)
    do ngcy = 1, miny
      if (ngcy.LE.nremy) then
        nrepy(ngcy) = noccy + 1
      else
        nrepy(ngcy) = noccy
      endif
    enddo
    ! Z-direction.
    noccz = set%ncoverz / groups%ngcellz
    nremz = set%ncoverz - noccz*groups%ngcellz
    minz  = min(set%ncoverz,groups%ngcellz)
    if (minz.GT.groups%mx_gedge) &
      call cq_abort('allocate_CS: too many groups in z-edgs', minz)
    do ngcz = 1, minz
      if (ngcz.LE.nremz) then
        nrepz(ngcz) = noccz + 1
      else
        nrepz(ngcz) = noccz
      endif
    enddo

    ! Go over groups in GCS periodic-irreducible set, calculating
    ! sim-cell (node-order, home-start) label of each.
    ng_in_cell = groups%ngcellx*groups%ngcelly*groups%ngcellz !No. of groups in sim-cell.
    ng_in_min  = minx*miny*minz
    ind = 0
    do ngcx = 1, minx
      do ngcy = 1, miny
        do ngcz = 1, minz
          ind = ind + 1
          ! x,y,z numbering of groups in sim-cell.
          nqx=1+mod(set%nx_origin+ngcx-set%nspanlx-2+nmodx, groups%ngcellx)
          nqy=1+mod(set%ny_origin+ngcy-set%nspanly-2+nmody, groups%ngcelly)
          nqz=1+mod(set%nz_origin+ngcz-set%nspanlz-2+nmodz, groups%ngcellz)
          ! CC in sim-cell.
          ind_qart=(nqx-1)*groups%ngcelly*groups%ngcellz+(nqy-1)*groups%ngcellz+nqz
          ! Node qwning ind_qart.
          ino=groups%inv_ngnode(ind_qart)
          ind_min(ind)=1+mod(ino-groups%inode_beg(inode)+ng_in_cell, ng_in_cell)
          ngcx_min(ind)=ngcx
          ngcy_min(ind)=ngcy
          ngcz_min(ind)=ngcz
        enddo
      enddo
    enddo
    ! sort minimum CS by nodes.
    call indexx(groups%mx_gcell,ng_in_min,ind_min,min_sort)

    ! Go over all GCS groups in NOPG order.
    ind_cover = 0
    do ind = 1, ng_in_min
      ngcx = ngcx_min(min_sort(ind))
      ngcy = ngcy_min(min_sort(ind))
      ngcz = ngcz_min(min_sort(ind))
      do nrx = 1, nrepx(ngcx)
        do nry = 1, nrepy(ngcy)
          do nrz = 1, nrepz(ngcz)
            ind_cover = ind_cover + 1
            nx_in_cover(ind_cover) = ngcx-1-set%nspanlx+(nrx-1)*groups%ngcellx
            ny_in_cover(ind_cover) = ngcy-1-set%nspanly+(nry-1)*groups%ngcelly
            nz_in_cover(ind_cover) = ngcz-1-set%nspanlz+(nrz-1)*groups%ngcellz
          enddo
        enddo
      enddo
    enddo

    ! Make CS for PS.
    nm_in_cover = 0
    do ind_cover = 1, set%ng_cover
      nsx = nx_in_cover(ind_cover)
      nsy = ny_in_cover(ind_cover)
      nsz = nz_in_cover(ind_cover)
      nqx = 1+mod(set%nx_origin+nsx+nmodx-1, groups%ngcellx)
      nqy = 1+mod(set%ny_origin+nsy+nmody-1, groups%ngcelly)
      nqz = 1+mod(set%nz_origin+nsz+nmodz-1, groups%ngcellz)
      ind_qart=(nqx-1)*groups%ngcelly*groups%ngcellz+(nqy-1)*groups%ngcellz+nqz ! CC in sim-cell.
      nm_in_cover = nm_in_cover + groups%nm_group(ind_qart)
    enddo
    set%mx_mcover = nm_in_cover

    ! Allocate arrays in set related to members.
    allocate (set%xcover(set%mx_mcover),set%ycover(set%mx_mcover),set%zcover(set%mx_mcover), &
              set%ig_cover(set%mx_mcover), STAT=stat_alloc)
    if (stat_alloc.NE.0) &
      call cq_abort('Error allocating arrays related to set:', set%mx_mcover)

    ! Deallocation.
    deallocate (nrepx,nrepy,nrepz, STAT=stat_alloc)
    if (stat_alloc.NE.0) &
      call cq_abort('Error deallocating nrepx,y,z:', groups%mx_gedge)
    deallocate (ind_min,ngcx_min,ngcy_min,ngcz_min,min_sort, STAT=stat_alloc)
    if (stat_alloc.NE.0) &
      call cq_abort('Error deallocating ind_min,ngcx,y,z_min,min_sort::', groups%mx_gcell)

    return
  end subroutine allocate_CSmember
  !!***

  !!****f* UpdateMember_module/cover_update_mparts
  !!
  !!  NAME 
  !!   cover_update_mparts
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Updates all the necessary arrays and variables in covering sets
  !!   for member-updates
  !!  INPUTS
  !!   set,groups,nx_in_cover,ny_in_cover,nz_in_cover,nmodx,nmody,nmodz,
  !!   dcellx,dcelly,dcellz,prim
  !!    - 'prim' is optional
  !!         |
  !!         |--> if 'set' is BCS_parts, 'prim' is 'bundle'
  !!         |--> if 'set' is DCS_parts, 'prim' is 'domain'
  !!         |--> if 'set' is 'ewald_CS' or 'D2_Cs', we don't need 'prim'
  !!  USES
  !!   global_module,basic_types
  !!  AUTHOR
  !!   Michiaki Arita
  !!  CREATION DATE
  !!   2013/07/02
  !!  MODIFICATION HISTORY
  !!
  !!  SOURCE
  !!
  subroutine cover_update_mparts(set,groups,nx_in_cover,ny_in_cover,nz_in_cover, &
                                 nmodx,nmody,nmodz,dcellx,dcelly,dcellz,prim)

    ! Module usage
    use basic_types
    use global_module, ONLY: rcellx,rcelly,rcellz,ni_in_cell, &
                             x_atom_cell,y_atom_cell,z_atom_cell, &
                             id_glob
    use GenComms, ONLY: cq_abort
    ! DB
    use io_module, ONLY: get_file_name
    use GenComms, ONLY: inode
    use global_module, ONLY: numprocs,io_lun
    use input_module, ONLY: io_assign,io_close

    implicit none

    ! passed variables
    integer :: nmodx,nmody,nmodz
    integer,allocatable :: nx_in_cover(:),ny_in_cover(:),nz_in_cover(:)
    real(double) :: dcellx,dcelly,dcellz
    type(cover_set)          :: set
    type(primary_set),target,optional :: prim
    type(group_set),target   :: groups

    ! local variables
    integer :: ind_cover,nsx,nsy,nsz,nqx,nqy,nqz,ind_qart
    integer :: ni,stat_alloc
    integer :: i,pr,ind_group,j
    real(double) :: xadd,yadd,zadd

    ! DB
    integer :: lun_db
    character(20) :: file_name

    ! Update n_ing_cover,icover_ibeg,xcover,ycover,zcover & ig_cover
    do ind_cover = 1, set%ng_cover              ! NOPG in CS.
      nsx = nx_in_cover(ind_cover)
      nsy = ny_in_cover(ind_cover)
      nsz = nz_in_cover(ind_cover)
      nqx = 1+mod(set%nx_origin+nsx+nmodx-1, groups%ngcellx)
      nqy = 1+mod(set%ny_origin+nsy+nmody-1, groups%ngcelly)
      nqz = 1+mod(set%nz_origin+nsz+nmodz-1, groups%ngcellz)
      xadd = (set%nx_origin+nsx-nqx)*dcellx
      yadd = (set%ny_origin+nsy-nqy)*dcelly
      zadd = (set%nz_origin+nsz-nqz)*dcellz
      ind_qart = set%lab_cell(ind_cover)        ! CC in sim-cell.
      if (groups%nm_group(ind_qart).GT.0) then
        set%n_ing_cover(ind_cover) = groups%nm_group(ind_qart)
        if (groups%icell_beg(ind_qart)+groups%nm_group(ind_qart)-1.GT.ni_in_cell) &
          call cq_abort('update_cover_mparts: x_atom_cover dim. exceeded', &
                        groups%icell_beg(ind_qart)+groups%nm_group(ind_qart)-1)
        do ni = 1, groups%nm_group(ind_qart)
          set%xcover(set%icover_ibeg(ind_cover)+ni-1) = &
            x_atom_cell(groups%icell_beg(ind_qart)+ni-1) + xadd
          set%ycover(set%icover_ibeg(ind_cover)+ni-1) = &
            y_atom_cell(groups%icell_beg(ind_qart)+ni-1) + yadd
          set%zcover(set%icover_ibeg(ind_cover)+ni-1) = &
            z_atom_cell(groups%icell_beg(ind_qart)+ni-1) + zadd
          set%ig_cover(set%icover_ibeg(ind_cover)+ni-1) = &
            id_glob(groups%icell_beg(ind_qart)+ni-1)
        enddo
        if (ind_cover.LT.set%ng_cover) &
          set%icover_ibeg(ind_cover+1)=set%icover_ibeg(ind_cover)+groups%nm_group(ind_qart)
      else ! If no atoms in ind_qart.
        if (ind_cover.LT.set%ng_cover) &
          set%icover_ibeg(ind_cover+1)=set%icover_ibeg(ind_cover)+groups%nm_group(ind_qart)
        set%n_ing_cover(ind_cover) = 0
      endif !(groups%nm_groups.GT.0)
    enddo

    ! Upate iprim_group.
    if (present(prim)) then
      pr = 0
      do i = 1, prim%groups_on_node
        nsx = prim%idisp_primx(i)
        nsy = prim%idisp_primy(i)
        nsz = prim%idisp_primz(i)
        ind_group = (nsx+set%nspanlx)*set%ncovery*set%ncoverz + &
                    (nsy+set%nspanly)*set%ncoverz + (nsz+set%nspanlz) + 1
        if (prim%nm_nodgroup(i).GT.0) then
          do j = 1, prim%nm_nodgroup(i)
            pr = pr + 1
            if (pr.GT.prim%mx_iprim) call cq_abort('cover_update_mparts: &
                                                    error in pr index', pr)
            set%iprim_group(pr) = ind_group
          enddo
        endif
      enddo
    endif

    ! Deallocation.
    deallocate (nx_in_cover,ny_in_cover,nz_in_cover, STAT=stat_alloc)
    if (stat_alloc.NE.0) &
      call cq_abort('Error deallocating nx_in_cover:', set%ng_cover)

    return
  end subroutine cover_update_mparts
  !!***

  !!****f* UpdateMember_module/updateMembers *
  !!
  !!  NAME 
  !!   updateMembers
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   All member-update processes are carried out
  !!  INPUTS
  !!   fixed_potential,velocity,iteration
  !!    - 'iteration' is optional
  !!    - 'iteration' will be deleted in the next update
  !!  USES
  !!   basic_types,global_module,ewald_module,group_module,primary_module,cover_module
  !!  AUTHOR
  !!   Michiaki Arita
  !!  CREATION DATE
  !!   2013/07/02
  !!  MODIFICATION HISTORY
  !!   2013/08/21 M.Arita
  !!   - Removed iteration
  !!  SOURCE
  !!
  !ORI subroutine updateMembers(fixed_potential,velocity,iteration)
  subroutine updateMembers(fixed_potential,velocity)

    ! Module usage
    use basic_types
    use global_module, ONLY: flag_dft_d2,ni_in_cell
    use GenComms, ONLY: cq_abort,my_barrier
    use ewald_module, ONLY: flag_old_ewald
    use group_module, ONLY: parts
    use primary_module, ONLY: bundle,domain
    use cover_module, ONLY: BCS_parts,DCS_parts,ewald_CS,D2_CS
    use density_module, ONLY: build_Becke_weights
    ! DB
    use global_module, ONLY: io_lun
    use input_module, ONLY: inode,ionode,io_assign,io_close
    use io_module, ONLY: get_file_name
    use global_module, ONLY: numprocs

    implicit none

    ! passed variables
    logical :: fixed_potential
    real(double) :: velocity(3,ni_in_cell)
    !ORI integer, intent(in), optional :: iteration

    ! local variables
    integer :: nmodx,nmody,nmodz
    integer, allocatable :: nx_in_cover(:),ny_in_cover(:),nz_in_cover(:)
    real(double) :: dcellx,dcelly,dcellz

    ! DB
    integer :: lun_db,lun_db2,lun_db3,lun_db4,lun_db5,lun_db6,lun_db7,lun_db8,lun_db9,lun_db10
    character(20) :: file_name

    ! NOTE: This if statement will be deleted in the next update
    !ORI if (present(iteration)) then                                ! for MD
    !ORI   call group_update_mparts(velocity,iteration)
    !ORI else
    !ORI   call group_update_mparts(velocity)                        ! for CG
    !ORI endif
    call group_update_mparts(velocity) ! both for md & cg

    ! Update members n PS
    call deallocate_PSmember
    call allocate_PSmember(parts)
    call primary_update_mparts
    call deallocate_CSmember(BCS_parts)
    call allocate_CSmember(BCS_parts,parts,nx_in_cover,ny_in_cover,nz_in_cover, &
                           nmodx,nmody,nmodz,dcellx,dcelly,dcellz,bundle)
    call cover_update_mparts(BCS_parts,parts,nx_in_cover,ny_in_cover,nz_in_cover, &
                             nmodx,nmody,nmodz,dcellx,dcelly,dcellz,bundle)
    if (inode.EQ.ionode) write (io_lun,*) "Finished updating BCS_parts!"
    ! Update members in CS
    call deallocate_CSmember(DCS_parts)
    call allocate_CSmember(DCS_parts,parts,nx_in_cover,ny_in_cover,nz_in_cover, &
                           nmodx,nmody,nmodz,dcellx,dcelly,dcellz,domain)
    call cover_update_mparts(DCS_parts,parts,nx_in_cover,ny_in_cover,nz_in_cover, &
                             nmodx,nmody,nmodz,dcellx,dcelly,dcellz,domain)
    if (inode.EQ.ionode) write (io_lun,*) "Finished updating DCS_parts!"
    if (.NOT. flag_old_ewald) then
      call deallocate_CSmember(ewald_CS)
      call allocate_CSmember(ewald_CS,parts,nx_in_cover,ny_in_cover,nz_in_cover, &
                             nmodx,nmody,nmodz,dcellx,dcelly,dcellz)
      call cover_update_mparts(ewald_CS,parts,nx_in_cover,ny_in_cover,nz_in_cover, &
                               nmodx,nmody,nmodz,dcellx,dcelly,dcellz)
      if (inode.EQ.ionode) write (io_lun,*) "Finished updating ewald_CS!"
    endif
    if (flag_dft_d2) then
      call deallocate_CSmember(D2_CS)
      call allocate_CSmember(D2_CS,parts,nx_in_cover,ny_in_cover,nz_in_cover, &
                             nmodx,nmody,nmodz,dcellx,dcelly,dcellz)
      call cover_update_mparts(D2_CS,parts,nx_in_cover,ny_in_cover,nz_in_cover, &
                               nmodx,nmody,nmodz,dcellx,dcelly,dcellz)
      if (inode.EQ.ionode) write (io_lun,*) "Finished updating D2_CS!"
    endif
    call my_barrier()
    !if (flag_Becke_weights) call build_Becke_weights

    return
  end subroutine updateMembers
  !! ***

end module UpdateMember_module
