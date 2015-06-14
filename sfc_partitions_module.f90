! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
!----------------------------------------------------------------------
! $Id$
!----------------------------------------------------------------------
! Module sfc_partitions_module
!----------------------------------------------------------------------
! Code area 1: Initialisation
!----------------------------------------------------------------------

!****h* Conquest/sfc_partitions_module
! PURPOSE
!   Module for creating partitions using Hilbert curve
! AUTHOR
!   L.Tong
! CREATION DATE
!   2013/04/03
! MODIFICATION HISTORY
! SOURCE
!
module sfc_partitions_module

  use datatypes

  implicit none

  save
  private

  public :: &
       sfc_partitions_to_processors, &
       n_parts_user, &
       average_atomic_diameter

  ! RCS tag for object file identification
  character(len=80), private :: &
       RCSid = "$Id$"

  ! partition type
  type part
     integer :: i_cc      ! cc index
     integer :: i_hilbert ! Hilbert curve index
     integer :: n_atoms   ! number of atoms it contains
     integer :: n_sf      ! number of support functions it contains
  end type part
  ! for recording size
  integer, parameter :: size_part = 4  ! in units of sizeof(integer)

  ! process type
  type proc
     integer :: id           ! processor id (rank within MPI communicator)
     integer :: i_part_start ! hilbert index of first partition on proc
     integer :: n_parts      ! number of partitions belonging to it
     integer :: n_atoms      ! number of atoms belonging to it
     integer :: n_sf         ! number of support functions belonging to it
  end type proc
  ! for recording size
  integer, parameter :: size_proc = 5 ! in units of sizeof(integer)

  ! cell type
  type cell
     real(double), dimension(3)   :: dims        ! cell dimensions
     real(double), dimension(3)   :: r_atoms_min ! minimum atomic coordinates
     real(double), dimension(3)   :: r_atoms_max ! maximum atomic coordinates
     real(double), dimension(2,3) :: gap         ! gap(1,i) is the start coordinate
                                                 ! of gap in i-th direction,
                                                 ! and gap(2,i) is the end coordinate
                                                 ! of gap in i-th direction
     integer               :: system_type    ! 3 = bulk, 2 = slab, 1 = chain, 0 = molecule
     logical, dimension(3) :: has_pbc        ! whether to treat direction to have pbc
     logical, dimension(3) :: is_folded      ! whether the struct is folded under pbc
     logical               :: is_translated  ! whether we have translated the original structure
     real(double), dimension(3) :: trans_vec ! translation vector of the
                                             ! atoms if the atoms were translated
  end type cell

  ! module variables
  real(double)               :: average_atomic_diameter ! used for determining gaps
  integer,      dimension(3) :: n_parts_user            ! user requested number of partitions
  real(double), parameter    :: no_pbc = 0.5_double     ! gap/cell_dim criteria
                                                     ! for no assuming pbc in
                                                     ! direction
  integer,      dimension(3) :: dim_nparts_fixed ! dims which n_parts is fixed
  integer,      dimension(3) :: dim_nparts_auto  ! dims which n_parts is not fixed
  integer :: n_dim_fixed      ! number of dims fixed
  integer :: n_dim_auto       ! numbet of dims not fixed
  integer,      dimension(3) :: n_parts ! number of partitions actually used
  integer,      dimension(3) :: n_divs  ! n_parts(i) = 2**n_divs(i)
  real(double), dimension(3) :: r_part  ! size of partition cells
  type(part),   dimension(:), allocatable :: parts ! information about partitions
  type(proc),   dimension(:), allocatable :: procs ! information about processors
  type(cell) :: FSC             ! information about the cell
  integer    :: max_natoms_part ! maximum allowed atoms per partition
  integer,    dimension(:,:), allocatable :: iglob_atom_part ! iglob_atom_part(i,j)
                                                             ! is global atomic
                                                             ! index of ith atom
                                                             ! in jth partition
!*****

contains

  !****f* sfc_partitions_module/sfc_partitions_to_processors
  ! PURPOSE
  !   From atomic structure, simulation cell shape, number of
  !   processors and user input, works out optimal number of
  !   partitions to devide the cell, and mapped to an 1D index using a
  !   flexible (non-cubic) 3D Hilbert curve. The partitions are then
  !   devided automatically among the processors optimising load
  !   balancing.
  ! USAGE
  !   call sfc_partitions_to_processors(cq_parts)
  ! INPUTS
  !   type(group_set) cqParts: the conquest partition group_set data
  ! OUTPUT
  !   setups and updates cqParts, effectively assigning computed
  !   partitions and atoms to processors
  ! AUTHOR
  !   L.Tong
  ! CREATION DATE
  !   2013/04/03
  ! MODIFICATION HISTORY
  !   2014/09/15 18:30 lat
  !    Fixed call start/stop_timer to timer_module (not timer_stdlocks_module !)
  !   2015/06/08 lat
  !    Added experimental backtrace
  ! SOURCE
  !
  subroutine sfc_partitions_to_processors(cqParts)

    use datatypes
    use global_module,          only: ni_in_cell, numprocs, species_glob, &
                                      iprint_init, io_lun, area_init
    use GenComms,               only: my_barrier, inode, ionode, cq_abort
    use species_module,         only: nsf_species
    use basic_types,            only: group_set
    use memory_module,          only: reg_alloc_mem, reg_dealloc_mem, type_int
    use timer_module,           only: start_timer,     stop_timer,    cq_timer
    use timer_module,           only: start_backtrace, stop_backtrace
    use timer_stdclocks_module, only: tmr_std_initialisation

    implicit none

    ! passed parameters
    type(group_set), intent(inout) :: cqParts ! CQ partition information type

    ! local variables
    logical :: reassign
    integer :: tot_natoms_parts, tot_n_sf, tot_n_parts
    integer :: ii, i_part, stat, counter
    integer :: min_nparts_proc, min_natoms_proc, min_nsf_proc
    integer ::  av_nparts_proc,  av_natoms_proc,  av_nsf_proc
    type(cq_timer) :: backtrace_timer

!****lat<$
    call start_backtrace(t=backtrace_timer,who='sfc_partitions_to_proc', &
                     where=1,level=3)
!****lat>$

    if ((inode == ionode) .and. (iprint_init > 2)) then
       write (io_lun, "(/,10x,a,/)") &
            "Using Hilbert curve automatic paritions"
    end if
    call start_timer(tmr_std_initialisation)
    ! allocate memory
    allocate(procs(numprocs), STAT=stat)
    if (stat /= 0) then 
       call cq_abort("sfc_partitions_to_processors(): Error alloating procs",&
                      numprocs)
    end if
    call reg_alloc_mem(area_init, numprocs * size_proc, type_int)
    ! setup partitions, this allocates parts
    call setup_partitions()
    ! check if the partitions are setup correctly
    call check_partitions()
    tot_n_parts = product(n_parts)
    tot_n_sf = 0
    do ii = 1, tot_n_parts
       tot_n_sf = tot_n_sf + parts(ii)%n_sf
    end do
    ! assign partitions to processors
    ! initial guess on the minimum amount of data to assign to each processor
    min_nparts_proc = tot_n_parts / numprocs
    min_natoms_proc = ni_in_cell / numprocs
    min_nsf_proc = tot_n_sf / numprocs
    ! the above must all >= 1, no need to check here however, as this
    ! must be the case if the partitions passed check_partitions()
    ! test.
    ! do load balancing
    reassign = .true.
    loop_reassign: do while (reassign)
       reassign = .false.
       ! loop over the processors, and assign partitions to each
       ! processor until nsf has reached or just past min_nsf_proc.
       ! dump all reminging partitions to last processor
       i_part = 1
       do ii = 1, numprocs
          ! initialise numbers
          procs(ii)%id = ii - 1
          procs(ii)%i_part_start = 0
          procs(ii)%n_parts = 0
          procs(ii)%n_atoms = 0
          procs(ii)%n_sf = 0
          counter = 1
          loop_part: do while (i_part <= tot_n_parts)
             if ((criteria(procs(ii)%n_parts, &
                           procs(ii)%n_atoms, &
                           procs(ii)%n_sf, &
                           min_nparts_proc, &
                           min_natoms_proc, &
                           min_nsf_proc)) &
                 .and. &
                 (ii < numprocs)) then
                exit loop_part
             end if
             if (counter == 1) procs(ii)%i_part_start = i_part
             procs(ii)%n_parts = procs(ii)%n_parts + 1
             procs(ii)%n_atoms = procs(ii)%n_atoms + parts(i_part)%n_atoms
             procs(ii)%n_sf = procs(ii)%n_sf + parts(i_part)%n_sf
             i_part = i_part + 1
             counter = counter + 1
          end do loop_part
       end do ! procs
       ! now that all processors has been assined partitions,
       ! check load balance, no need for single processor calculations
       if (numprocs > 1) then
          if (criteria(min_nparts_proc, &
                       min_natoms_proc, &
                       min_nsf_proc, &
                       procs(numprocs)%n_parts, &
                       procs(numprocs)%n_atoms, &
                       procs(numprocs)%n_sf)) then
             ! too few data allocated in the last processor
             reassign = .true.
             min_nparts_proc = min_nparts_proc - 1
             min_natoms_proc = min_natoms_proc - 1
             min_nsf_proc = min_nsf_proc - 1
             cycle loop_reassign
          end if
       end if
    end do loop_reassign
    ! the assignment loop cycles until there are too much data
    ! in the last processor, at which point, shift the partition
    ! staring points on some of the processors
    if (numprocs > 1) then
       av_nparts_proc = nint(real(tot_n_parts,double) / real(numprocs,double))
       av_natoms_proc = nint(real(ni_in_cell,double) / real(numprocs,double))
       av_nsf_proc = nint(real(tot_n_sf,double) / real(numprocs,double))
       do ii = numprocs, 2, -1
          do while (criteria((procs(ii)%n_parts - av_nparts_proc + 1), &
                             (procs(ii)%n_atoms - av_natoms_proc + 1), &
                             (procs(ii)%n_sf - av_nsf_proc + 1), &
                             1, &
                             parts(procs(ii)%i_part_start)%n_atoms, &
                             parts(procs(ii)%i_part_start)%n_sf))
             ! move one partition from proc ii to ii-1
             i_part = procs(ii)%i_part_start
             procs(ii)%i_part_start = procs(ii)%i_part_start + 1
             procs(ii)%n_parts = procs(ii)%n_parts - 1
             procs(ii)%n_atoms = procs(ii)%n_atoms - parts(i_part)%n_atoms
             procs(ii)%n_sf = procs(ii)%n_sf - parts(i_part)%n_sf
             procs(ii-1)%n_parts = procs(ii-1)%n_parts + 1
             procs(ii-1)%n_atoms = procs(ii-1)%n_atoms + parts(i_part)%n_atoms
             procs(ii-1)%n_sf = procs(ii-1)%n_sf + parts(i_part)%n_sf
          end do
       end do
    end if
    ! reshaffle the empty partitions
    if (numprocs > 1) then
       do ii = 2, numprocs
          do while ((parts(procs(ii)%i_part_start)%n_atoms == 0) .and. &
                    (procs(ii)%n_parts > av_nparts_proc))
             ! move one partition from proc ii to ii-1
             i_part = procs(ii)%i_part_start
             procs(ii)%i_part_start = procs(ii)%i_part_start + 1
             procs(ii)%n_parts = procs(ii)%n_parts - 1
             procs(ii-1)%n_parts = procs(ii-1)%n_parts + 1
          end do
       end do
    end if
    ! construct Conquest partition data
    call construct_cq_partitions(cqParts)
    ! print information
    call print_information()
    call my_barrier()
    ! deallocate all the arrays
    if (allocated(parts)) then
       deallocate(parts, STAT=stat)
       if (stat /= 0) then
          call cq_abort("sfc_partitions_to_processors(): &
                         &Error dealloating parts")
       end if
       call reg_dealloc_mem(area_init, numprocs * size_part, type_int)
    end if
    if (allocated(procs)) then
       deallocate(procs, STAT=stat)
       if (stat /= 0) then
          call cq_abort("sfc_partitions_to_processors(): &
                         &Error dealloating parts")
       end if
       call reg_dealloc_mem(area_init, tot_n_parts * size_proc, &
                            type_int)
    end if
    if (allocated(iglob_atom_part)) then
       deallocate(iglob_atom_part, STAT=stat)
       if (stat /= 0) then
          call cq_abort("sfc_partitions_to_processors(): &
                         &Error dealloating iglob_atom_part")
       end if
       call reg_dealloc_mem(area_init, max_natoms_part * tot_n_parts, &
                            type_int)
    end if
    call stop_timer(tmr_std_initialisation)

!****lat<$
    call stop_backtrace(t=backtrace_timer,who='sfc_partitions_to_processors')
!****lat>$

    return
  end subroutine sfc_partitions_to_processors
  !*****


  !****f* sfc_partitions_module/criteria
  ! PURPOSE
  !   does comparison between two sets of test data with, if test* >
  !   crit*, then returns true, otherwise false. The aim of this
  !   function is to provide flexibility on which pair of Parts, Atoms
  !   or Sf data to compare, depending on the setting on global_module
  !   variable load_balance
  ! USAGE
  !   if (criteria(testParts, testAtoms, testSf, &
  !       critParts, critAtoms, critSf) then
  !     ...
  !   end if
  ! INPUTS
  !   integer testParts, critParts: first pair for comparison
  !   integer testAtoms, critAtoms: second pair for comparison
  !   integer testSf, critSf: third pair for comparison
  ! RETURN VALUE
  !   depending on the value of load_balance, compares only one pair
  !   of the three pairs, and returns true if test* > crit*, otherwise
  !   returns false
  ! AUTHOR
  !   L.Tong
  ! CREATION DATE
  !   2013/04/03
  ! MODIFICATION HISTORY
  ! SOURCE
  !
  function criteria(testParts, testAtoms, testSf, &
                    critParts, critAtoms, critSf)

    use global_module, only: load_balance

    implicit none

    ! passed parameters
    integer, intent(in) :: testParts, testAtoms, testSf
    integer, intent(in) :: critParts, critAtoms, critSf

    ! returned value
    logical :: criteria
    criteria = .false.
    select case (load_balance)
    case (0) ! balance with respect to number of atoms in proc
       if (testAtoms > critAtoms) criteria = .true.
    case (1) ! balance with respect to number of parts in proc
       if (testParts > critParts) criteria = .true.
    case (2) ! balance with respect to number of sf in proc
       if (testSf > critSf) criteria = .true.
    end select
  end function criteria
  !*****


  !****f* sfc_partitions_module/construct_cq_partitions
  ! PURPOSE
  !   construct the Conquest group_set data for partitions, using the
  !   computed module variables from this module
  ! USAGE
  !   call construct_cq_partitions(cq_parts)
  ! INPUTS
  !   type(group_set) cqParts: the conquest group_set data for partitions
  ! OUTPUT
  !   cqParts updated
  ! AUTHOR
  !   L.Tong
  ! CREATION DATE
  !   2013/04/03
  ! MODIFICATION HISTORY
  ! SOURCE
  !
  subroutine construct_cq_partitions(cqParts)

    use global_module,    only: numprocs, atom_coord, ni_in_cell, &
                                id_glob, id_glob_inv
    use basic_types,      only: group_set
    use maxima_module,    only: maxpartsproc, maxatomspart, &
                                maxatomsproc, maxpartscell
    use construct_module, only: init_group

    implicit none

    ! passed parameters
    type(group_set), intent(inout) :: cqParts

    ! local variables
    integer, dimension(3)                :: i_part_xyz
    integer, dimension(:,:), allocatable :: i_atom_global
    integer, dimension(:),   allocatable :: i_atom_part
    integer :: ii, ia, maxpartsedge, i_cc, counter, stat

    maxpartscell = product(n_parts)
    maxpartsedge = maxval(n_parts)
    maxatomsproc = 0
    maxpartsproc = 0
    do ii = 1, numprocs
       maxatomsproc = max(maxatomsproc, procs(ii)%n_atoms)
       maxpartsproc = max(maxpartsproc, procs(ii)%n_parts)
    end do
    maxatomspart = 0
    do ii = 1, maxpartscell
       maxatomspart = max(maxatomspart, parts(ii)%n_atoms)
    end do
    ! setup group_set cqParts data type
    call init_group(cqParts, maxpartsproc, maxpartsedge, &
                    maxpartscell, maxatomspart, numprocs)
    cqParts%ngcellx = n_parts(1)
    cqParts%ngcelly = n_parts(2)
    cqParts%ngcellz = n_parts(3)
    do ii = 1, numprocs
       cqParts%ng_on_node(ii) = procs(ii)%n_parts
       cqParts%inode_beg(ii) = procs(ii)%i_part_start
    end do
    counter = 1
    do ii = 1, maxpartscell
       i_cc = parts(ii)%i_cc
       cqParts%ngnode(ii) = i_cc
       cqParts%nm_group(i_cc) = parts(ii)%n_atoms
       cqParts%icell_beg(i_cc) = counter
       counter = counter + parts(ii)%n_atoms
       if (parts(ii)%n_atoms > 0) then
          do ia = 1, parts(ii)%n_atoms
             id_glob(cqParts%icell_beg(i_cc)+ia-1) = iglob_atom_part(ia,ii)
             id_glob_inv(iglob_atom_part(ia,ii)) = cqParts%icell_beg(i_cc)+ia-1
          end do
       end if
    end do

    return
  end subroutine construct_cq_partitions
  !*****


  !****f* sfc_partitions_module/setup_partitions
  ! PURPOSE
  !   Works out cell information, and setup partitions accordingly,
  !   then assign atoms to the partitions, and refine partitions if
  !   necessary to get ideal number of atoms in each partition
  ! USAGE
  !   call setup_partitions()
  ! INPUTS
  !   module variables:
  ! OUTPUT
  !   updates module variables
  ! AUTHOR
  !   L.Tong
  ! CREATION DATE
  !   2013/04/02
  ! MODIFICATION HISTORY
  !   2013/08/25 M.Arita
  !   - Introduced shift_in_bohr to make partitioning consistent to
  !     shift in matrix reconstruction
  ! SOURCE
  !
  subroutine setup_partitions()
    ! given nParts and cell, workout partition indices and number of
    ! atoms in partitions etc.
    ! parts should have dimension (2**nDivs(1) * 2**nDivs(2) * 2**nDivs(3))
    ! atoms should have dimension number of atoms in FSC
    ! parts array are assumed to be indexed by hilbert indices
    use datatypes
    use numbers
    use GenComms,       only: cq_abort
    use global_module,  only: atom_coord, ni_in_cell, species_glob, &
                              area_init, shift_in_bohr
    use species_module, only: nsf_species
    use Hilbert3D,      only: Hilbert3D_Initialise, Hilbert3D_IntTOCoords, &
                              Hilbert3D_CoordsToInt
    use memory_module,  only: reg_alloc_mem, reg_dealloc_mem, type_int

    implicit none

    ! local variables
    integer                    :: ii, stat
    integer                    :: ihilbert, icc, ia
    real(double)               :: v_part, area
    real(double), dimension(3) :: r_occupied_cell
    integer,      dimension(3) :: ipart_xyz, nParts
    integer                    :: n_parts_total
    logical                    :: refine

    

    ! get cell and system information
    call get_cell_info()
    ! get initial guess for number of partitions
    call get_initial_partitions()
    ! assign atoms to the partitions and find the optimal number of
    ! partitions based on the max number of atoms in each partition
    refine = .true.
    loop_refine: do while (refine)
       refine = .false.
       ! allocate  array
       n_parts_total = product(n_parts)
       allocate(parts(n_parts_total), STAT=stat)
       if (stat /= 0) then
          call cq_abort("setup_partitions: Error allocating parts", &
                        n_parts_total)
       end if
       call reg_alloc_mem(area_init, n_parts_total, type_int)
       allocate(iglob_atom_part(max_natoms_part,n_parts_total), STAT=stat)
       if (stat /= 0) then
          call cq_abort("setup_partitions: Error allocating parts", &
                        max_natoms_part, n_parts_total)
       end if
       call reg_alloc_mem(area_init, max_natoms_part * n_parts_total, &
                          type_int)
       ! initialise number of atoms and support functions
       iglob_atom_part = 0
       do ihilbert = 1, n_parts_total
          parts(ihilbert)%n_atoms = 0
          parts(ihilbert)%n_sf = 0
       end do
       call Hilbert3D_Initialise(n_divs(1), n_divs(2), n_divs(3))
       ! get dimensions of each partition cell
       r_part(1:3) = FSC%dims(1:3) / real(n_parts(1:3),double)
       ! loop over atoms in cell, work out which partitions they are
       ! in, and do the counting
       loop_atoms: do ia = 1, ni_in_cell
          ! get the partition x,y,z indices for the partition containing the atoms
          !ORI ipart_xyz(1:3) = floor(atom_coord(1:3,ia) / r_part(1:3) + RD_ERR)
          ipart_xyz(1:3) = floor((atom_coord(1:3,ia)+shift_in_bohr) / r_part(1:3) )
          ! get the corresponding Hilbert curve index corresponding to
          ! the partition
          call Hilbert3D_CoordsToInt(ipart_xyz, ihilbert)
          ! accumulate the atoms count
          parts(ihilbert)%n_atoms = parts(ihilbert)%n_atoms + 1
          ! check if the number of atoms per partition exceeds number allowed
          ! and if the partitioning needs to be refined
          if ((parts(ihilbert)%n_atoms > max_natoms_part) .and. &
              (n_dim_auto > 0) .and. &
              (minval(r_part) > 1.0_double)) then
             refine = .true.
             ii = refine_direction(r_part, n_dim_auto, dim_nparts_auto)
             n_divs(ii) = n_divs(ii) + 1
             n_parts(ii) = 2**n_divs(ii)
             r_part(ii) = FSC%dims(ii) / real(n_parts(ii),double)
             deallocate(parts, STAT=stat)
             if (stat /= 0) then
                call cq_abort("setup_partitions: Error deallocating parts")
             end if
             call reg_dealloc_mem(area_init, n_parts_total, type_int)
             deallocate(iglob_atom_part, STAT=stat)
             if (stat /= 0) then
                call cq_abort("setup_partitions: Error deallocating parts")
             end if
             call reg_dealloc_mem(area_init, max_natoms_part * n_parts_total, &
                                  type_int)
             cycle loop_refine
          end if
          ! accumulate the number of support functions
          parts(ihilbert)%n_sf = parts(ihilbert)%n_sf + &
                                 nsf_species(species_glob(ia))
          ! record the atomic number
          iglob_atom_part(parts(ihilbert)%n_atoms,ihilbert) = ia
       end do loop_atoms
       ! after obtaining the optimal partitioning levels, get rest
       ! partition information
       do ihilbert = 1, n_parts_total
          call Hilbert3D_IntToCoords(ihilbert, ipart_xyz)
          ! get cc index (start counting from 1)
          icc = ipart_xyz(1) * n_parts(2) * n_parts(3) + &
                ipart_xyz(2) * n_parts(3) + &
                ipart_xyz(3) + 1
          ! set relevent parts data
          parts(ihilbert)%i_cc = icc
          parts(ihilbert)%i_hilbert = ihilbert
       end do ! ihilbert
    end do loop_refine
    return
  end subroutine setup_partitions
  !*****


  !****f* sfc_partitions_module/get_initial_partitions
  ! PURPOSE
  !   Make a guess on the initial number of partitions based on user
  !   input and cell information obtained from get_cell_info
  ! USAGE
  !   call get_initial_partitions()
  ! INPUTS
  !   reads module variables:
  !   FSC, global_maxatomspart, n_parts_user, atom_coord, numprocs,
  !   ni_in_cell
  ! OUTPUT
  !   updates module variables:
  !   n_parts, r_part, dim_nparts_fixed, dim_nparts_auto, n_dim_fixed,
  !   n_dim_auto, n_divs, max_natoms_part
  ! AUTHOR
  !   L.Tong
  ! CREATION DATE
  !   2013/04/02
  ! MODIFICATION HISTORY
  ! SOURCE
  !
  subroutine get_initial_partitions()

    use datatypes
    use numbers
    use global_module, only: atom_coord, global_maxatomspart, &
                             ni_in_cell, numprocs
    implicit none

    ! local variables
    real(double), dimension(3) :: r_occupied_cell
    real(double)               :: v_part, area
    integer                    :: ii
    ! get an estimate of maximum number of atoms per partition,
    ! minimum one partition per processor, ideally each partition
    ! should contain no more than ni_in_cell / numprocs number of
    ! atoms for good load balancing between processors
    max_natoms_part = ni_in_cell / numprocs
    ! should not below 1 or exceed the global limit set by user
    max_natoms_part = max(max_natoms_part, 1)
    max_natoms_part = min(max_natoms_part, global_maxatomspart)
    ! determine if the number of partitions in particular direction is
    ! fixed or to be determined automatically
    n_dim_fixed = 0
    n_dim_auto = 0
    dim_nparts_fixed = 0
    dim_nparts_auto = 0
    do ii = 1, 3
       if (n_parts_user(ii) == 0) then
          ! automatically determine number of partitions
          if (FSC%has_pbc(ii) .or. (FSC%system_type == 0)) then
             n_dim_auto = n_dim_auto + 1
             dim_nparts_auto(n_dim_auto) = ii
          else
             n_dim_fixed = n_dim_fixed + 1
             dim_nparts_fixed(n_dim_fixed) = ii
             n_parts(ii) = 1
             n_divs(ii) = 0
             r_part(ii) = FSC%dims(ii)
          end if
       else
          ! number of partitions set by the user
          n_dim_fixed = n_dim_fixed + 1
          dim_nparts_fixed(n_dim_fixed) = ii
          n_divs(ii) = int_log2(n_parts_user(ii))
          n_parts(ii) = 2**n_divs(ii)
          r_part(ii) = FSC%dims(ii) / real(n_parts(ii),double)
       end if
    end do
    ! if automatic determination of number of partitions is required
    if (n_dim_auto > 0) then
       do ii = 1, 3
          if (FSC%is_folded(ii)) then
             r_occupied_cell(ii) = FSC%dims(ii) - FSC%gap(2,ii) + FSC%gap(1,ii)
          else
             r_occupied_cell(ii) = FSC%r_atoms_max(ii) - FSC%r_atoms_min(ii)
          end if
       end do
       v_part = (real(max_natoms_part,double) / real(ni_in_cell,double)) * &
                product(r_occupied_cell)
       select case (n_dim_auto)
       case (3)
          do ii = 1, 3
             r_part(ii) = v_part**third
          end do
       case (2)
          area = v_part / &
                 min(r_part(dim_nparts_fixed(1)), &
                     r_occupied_cell(dim_nparts_fixed(1)))
          do ii = 1, 2
             r_part(dim_nparts_auto(ii)) = sqrt(area)
          end do
       case (1)
          r_part(dim_nparts_auto(1)) = &
               v_part / (min(r_part(dim_nparts_fixed(1)), &
                             r_occupied_cell(dim_nparts_fixed(1))) * &
                         min(r_part(dim_nparts_fixed(2)), &
                             r_occupied_cell(dim_nparts_fixed(2))))
       end select
       do ii = 1, 3
          n_parts(ii) = max(nint(FSC%dims(ii) / r_part(ii)), 1)
          n_divs(ii) = int_log2(n_parts(ii))
          n_parts(ii) = 2**n_divs(ii)
          r_part(ii) = FSC%dims(ii) / real(n_parts(ii),double)
       end do
    end if

    return
  end subroutine get_initial_partitions
  !*****


  !****f* sfc_partitions_module/int_log2
  ! PURPOSE
  !   returns the minimum y, where 2**y >= x
  ! USAGE
  !   y = int_log2(x)
  ! INPUTS
  !   integer x
  ! RETURN VALUE
  !   integer y
  ! AUTHOR
  !   L.Tong
  ! CREATION DATE
  !   2013/03/30
  ! MODIFICATION HISTORY
  ! SOURCE
  !
  function int_log2(x)
    ! returns the minumum power (of 2) that is greater or equal to num

    implicit none

    ! passed parameters
    integer, intent(in) :: x

    ! returned value
    integer :: int_log2
    int_log2 = 0
    do while (2**int_log2 < x)
       int_log2 = int_log2 + 1
    end do

    return
  end function int_log2
  !*****


  !****f* sfc_partitions_module/
  ! PURPOSE
  !   Works out the direction along which to refine the partition divisions
  !   Rule:
  !   1) If the number of partitions are fixed, do not refine
  !   2) Refine the direction which has the longest partition cell length
  ! USAGE
  !   dim = refine_direction(rParts, fixedParts)
  ! INPUTS
  !   integer rParts(ndims) : length in idim-th direction of each partition
  !   integer nRefineDims : number of dimensions along which partitioning
  !                         can be refined
  !   integer refineDims(3) : the first nRefineDims elements stores the
  !                           dimensions along which partitioning can be
  !                           refined
  ! RETURN VALUE
  !   integer refine_direction : direction along which to refine
  ! AUTHOR
  !   L.Tong
  ! CREATION DATE
  !   2013/03/30
  ! MODIFICATION HISTORY
  ! SOURCE
  !
  function refine_direction(rParts, nRefineDims, refineDims)

    use datatypes

    implicit none

    ! passed parameters
    real(double), dimension(:), intent(in) :: rParts
    integer,                    intent(in) :: nRefineDims
    integer,      dimension(:), intent(in) :: refineDims

    ! returned value
    integer :: refine_direction

    ! local variables
    integer      :: ii, dim
    real(double) :: max_r

    max_r = 0.0_double
    do ii = 1, nRefineDims
       dim = refineDims(ii)
       if (rParts(dim) > max_r) then
          max_r = rParts(dim)
          refine_direction = dim
       end if
    end do

  end function refine_direction
  !*****


  !****f* sfc_partitions_module/get_cell_info
  ! PURPOSE
  !   Work out the type of system (bulk, slab, chain or molecule), the
  !   extent of the atoms, and if there are any gaps exist in the
  !   structure
  ! USAGE
  !   call get_cell_info
  ! INPUTS
  !   Data from module variables and atomic coordinates
  ! OUTPUT
  !   FSC module variable updated
  ! AUTHOR
  !   L.Tong
  ! CREATION DATE
  !   2013/04/01
  ! MODIFICATION HISTORY
  !   2013/08/25 M.Arita
  !   - Introduced shift_in_bohr and made the code consistent to shift in
  !     matrix reconstruction
  ! SOURCE
  !
  subroutine get_cell_info()

    use datatypes
    use numbers
    use dimens,        only: r_super_x, r_super_y, r_super_z
    use global_module, only: atom_coord, ni_in_cell, area_init, shift_in_bohr
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_int
    use GenComms,      only: cq_abort

    implicit none

    ! local variables
    integer,      parameter    :: min_n_slices = 10
    integer,      dimension(3) :: n_blocks
    integer,      dimension(3) :: i_block_xyz
    real(double), dimension(3) :: r_block
    real(double), dimension(2) :: limits
    integer,      dimension(:), allocatable :: n_atoms_block
    integer      :: iatom, n_blocks_total, dim_min_r_block
    integer      :: ii, icc, stat
    integer      :: gap_block_min, gap_block_max
    real(double) :: r_gap_block_min, r_gap_block_max
    real(double) :: tmp_min, tmp_max
    real(double) :: min_r_block

    ! get cell dimensions
    FSC%dims(1) = r_super_x
    FSC%dims(2) = r_super_y
    FSC%dims(3) = r_super_z
    ! find the mimimum r_block
    r_block(1) = r_super_x / real(min_n_slices, double)
    r_block(2) = r_super_y / real(min_n_slices, double)
    r_block(3) = r_super_z / real(min_n_slices, double)
    min_r_block = r_block(1)
    dim_min_r_block = 1
    do ii = 2, 3
       if (r_block(ii) < min_r_block) then
          min_r_block = r_block(ii)
          dim_min_r_block = ii
       end if
    end do
    ! the minimum r_block should not be below average_atomic_diameter
    if (min_r_block < average_atomic_diameter) then
       min_r_block = average_atomic_diameter
    end if
    ! recalculate the r_blocks based on the minimum value
    n_blocks(1) = nint(r_super_x / min_r_block)
    n_blocks(2) = nint(r_super_y / min_r_block)
    n_blocks(3) = nint(r_super_z / min_r_block)
    r_block(1) = r_super_x / real(n_blocks(1), double)
    r_block(2) = r_super_y / real(n_blocks(2), double)
    r_block(3) = r_super_z / real(n_blocks(3), double)
    n_blocks_total = n_blocks(1) * n_blocks(2) * n_blocks(3)
    allocate(n_atoms_block(n_blocks_total), STAT=stat)
    if (stat /= 0) then
       call cq_abort("get_cell_info: Error allocating n_atoms_block", &
                     n_blocks_total)
    end if
    call reg_alloc_mem(area_init, n_blocks_total, type_int)
    n_atoms_block = 0.0_double
    ! calculate the number of atoms in each block
    do iatom = 1, ni_in_cell
       !ORI i_block_xyz(1:3) = floor(atom_coord(1:3,iatom) / r_block(1:3) + RD_ERR)
       i_block_xyz(1:3) = floor((atom_coord(1:3,iatom)+shift_in_bohr) / r_block(1:3) ) 
       icc = i_block_xyz(1) * n_blocks(2) * n_blocks(3) + &
             i_block_xyz(2) * n_blocks(3) + &
             i_block_xyz(3) + 1
       n_atoms_block(icc) = n_atoms_block(icc) + 1
    end do
    ! check for gaps and determine system type
    FSC%system_type = 3
    FSC%has_pbc = .true.
    FSC%is_folded = .false.
    FSC%gap = 0.0_double
    do ii = 1, 3
       ! get extent of atoms
       call get_extent(ii, atom_coord, FSC%r_atoms_min(ii), &
                       FSC%r_atoms_max(ii))
       call check_gap(ii, n_blocks, n_atoms_block, &
                      gap_block_min, gap_block_max)
       ! if gaps are found
       if ((gap_block_min > 0) .and. (gap_block_max > 0)) then
          ! get gap extent, remember block_gap counts from 1
          r_gap_block_min = (gap_block_min - 1) * r_block(ii)
          r_gap_block_max = (gap_block_max - 1) * r_block(ii)
          if ((r_gap_block_min >= FSC%r_atoms_min(ii)) .and. &
               (r_gap_block_max <= FSC%r_atoms_max(ii))) then
             ! the gap is in the middle of the cell
             FSC%is_folded(ii) = .true.
             limits(1) = FSC%r_atoms_min(ii)
             limits(2) = r_gap_block_min
             call get_extent(ii, atom_coord, tmp_min, FSC%gap(1,ii), limits)
             limits(1) = r_gap_block_max
             limits(2) = FSC%r_atoms_max(ii)
             call get_extent(ii, atom_coord, FSC%gap(2,ii), tmp_max, limits)
             ! check if no pbc in this direction
             if ((FSC%gap(2,ii) - FSC%gap(1,ii)) / FSC%dims(ii) > no_pbc) then
                FSC%system_type = FSC%system_type - 1
                FSC%has_pbc(ii) = .false.
             end if
          else
             ! the gap is either on top or underneath the atoms
             ! check if no pbc in this direction
             if ((FSC%r_atoms_max(ii) - FSC%r_atoms_min(ii)) / FSC%dims(ii) <= &
                 no_pbc) then
                FSC%system_type = FSC%system_type - 1
                FSC%has_pbc(ii) = .false.
             end if
             if (r_gap_block_min >= FSC%r_atoms_max(ii)) then
                ! the gap is on top of the atoms
                FSC%gap(1,ii) = FSC%r_atoms_max(ii)
                FSC%gap(2,ii) = FSC%dims(ii)
             else if (r_gap_block_max <= FSC%r_atoms_min(ii)) then
                ! the gap is underneath the atoms
                FSC%gap(1,ii) = 0.0_double
                FSC%gap(2,ii) = FSC%r_atoms_min(ii)
             end if
          end if
       end if
    end do ! ii
    deallocate(n_atoms_block, STAT=stat)
    if (stat /= 0) then
       call cq_abort("get_cell_info: Error deallocating n_atoms_block")
    end if
    call reg_dealloc_mem(area_init, n_blocks_total, type_int)

    return
  end subroutine get_cell_info
  !*****


  !****f* sfc_partitions_module/get_extent
  ! PURPOSE
  !   Get the extent of a collection of atomic coordinates in a given
  !   direction, within an optional limit
  ! USAGE
  !   - To find the absolute extent:
  !     call get_extent(dir, atom_coords, min, max)
  !   - To find the extent within a given limit
  !     call get_extent(dir, atom_coords, min, max, limits)
  ! INPUTS
  !   integer dir : direction, dir = 1, 2 or 3
  !   real(double) atom_coords(3,natoms) : atomic coords
  !   real(double), optional limits(2) : find extent within [limits(1), limits(2)]
  ! OUTPUT
  !   real(double) extent_min : minimum coordinate in dir direction
  !   real(double) extent_max : maximum coordinate in dir direction
  ! AUTHOR
  !   L.Tong
  ! CREATION DATE
  !   2013/03/31
  ! MODIFICATION HISTORY
  ! SOURCE
  !
  subroutine get_extent(dir, atom_coords, extent_min, extent_max, limits)
    ! from atom_coords, get the minimum and maximum coordinates along
    ! direction dir, and within the given limits(1) and limits(2)

    use datatypes

    implicit none

    ! passed parameters
    integer,                       intent(in)          :: dir
    real(double), dimension(:,:), intent(in)           :: atom_coords
    real(double), dimension(2),   intent(in), optional :: limits
    real(double),                 intent(out)          :: extent_min, extent_max

    ! local variables
    integer :: ii
    integer :: n_atoms

    n_atoms    = size(atom_coords, 2)
    extent_min = huge(1.0_double)
    extent_max = 0.0_double
    if (present(limits)) then
       do ii = 1, n_atoms
          if ((atom_coords(dir,ii) > extent_max) .and. &
              (atom_coords(dir,ii) <= limits(2))) then
             extent_max = atom_coords(dir,ii)
          end if
          if ((atom_coords(dir,ii) < extent_min) .and. &
              (atom_coords(dir,ii) >= limits(1))) then
             extent_min = atom_coords(dir,ii)
          end if
       end do
    else
       do ii = 1, n_atoms
          if (atom_coords(dir,ii) > extent_max) then
             extent_max = atom_coords(dir,ii)
          end if
          if (atom_coords(dir,ii) < extent_min) then
             extent_min = atom_coords(dir,ii)
          end if
       end do
    end if
    ! in the case where extend is tiny, avoid rounding errors
    if (extent_min > extent_max)  extent_min = extent_max
    return
  end subroutine get_extent
  !*****


  !****f* sfc_partitions_module/check_gap
  ! PURPOSE
  !   Find the maximum gap along direction a. A gap corresponds to the
  !   slab of blocks normal to direction a that does not contain
  !   atoms. If gaps are found, the maximum gap is returned, with
  !   gap_start and gap_end containing the start index and end index
  !   of the blocks in direction a, for the largest empty slab. If no
  !   gaps are found, then gap_start and gap_end will be zero
  ! USAGE
  !   call check_gap(a, b, c, n_blocks, n_atoms_block, gap_start, gap_end)
  ! INPUTS
  !   integer dir : search direction for gap, dir = 1, 2, or 3
  !   integer n_blocks(3) : number of blocks in each of three directions
  !   integer n_atoms_block(tot_n_blocks) : number of atoms in each block
  ! OUTPUT
  !   integer gap_start : starting a index of empty slab (gap)
  !   integer gap_end   : end a index of empty slab (gap). The number of
  !                       blocks along a in the gap would be
  !                       gap_end - gap_start + 1
  ! AUTHOR
  !   L.Tong
  ! CREATION DATE
  !   2013/03/31
  ! MODIFICATION HISTORY
  ! SOURCE
  !
  subroutine check_gap(dir, n_blocks, n_atoms_block, gap_start, gap_end)
    ! check gap in direction dir, transverse directions are b and c
    ! find the largest gap

    implicit none

    ! passed parameters
    integer,               intent(in)  :: dir
    integer, dimension(:), intent(in)  :: n_blocks, n_atoms_block
    integer,               intent(out) :: gap_start, gap_end

    ! local variables
    integer, dimension(2) :: tran_dir ! transverse directions
    integer, dimension(3) :: cc_multiplier
    integer :: ii, ia, ib, ic, icc
    integer :: tmp_start, tmp_gap, gap
    logical :: start_count, found_atom

    ! work out transverse directions
    do ii = 1, 2
       tran_dir(ii) = mod(dir+ii-1, 3) + 1
    end do
    ! define cc index multipliers
    cc_multiplier(1) = n_blocks(2) * n_blocks(3)
    cc_multiplier(2) = n_blocks(3)
    cc_multiplier(3) = 1
    ! initialise values
    gap = 0
    start_count = .false.
    gap_start = 0
    gap_end = 0
    ! find the maximum gap
    loop_a: do ia = 1, n_blocks(dir)
       found_atom = .false.
       loop_b: do ib = 1, n_blocks(tran_dir(1))
          loop_c: do ic = 1, n_blocks(tran_dir(2))
             icc = (ia-1) * cc_multiplier(dir) + &
                   (ib-1) * cc_multiplier(tran_dir(1)) + &
                   (ic-1) * cc_multiplier(tran_dir(2)) + 1
             if (n_atoms_block(icc) > 0)  then
                found_atom = .true.
                ! we have reached the end of the gap
                if (start_count) then
                   start_count = .false.
                   tmp_gap = ia - tmp_start
                   if (tmp_gap > gap) then
                      gap = tmp_gap
                      gap_start = tmp_start
                      gap_end = ia - 1
                   end if
                end if
                exit loop_b
             end if
          end do loop_c
       end do loop_b
       ! only start counting after found an empty slab
       if ((.not. start_count) .and. (.not. found_atom)) then
          tmp_start = ia
          start_count = .true.
       end if
       ! if no atoms were found, but have reached the end of cell
       if ((start_count) .and. (.not. found_atom) .and. &
           (ia == n_blocks(dir))) then
          ! remember the slab is still empty, so must be counted in
          ! the gap, so need to add 1
          tmp_gap = ia - tmp_start + 1
          if (tmp_gap > gap) then
             gap = tmp_gap
             gap_start = tmp_start
             gap_end = ia
          end if
       end if
    end do loop_a

    return
  end subroutine check_gap
  !*****


  !****f* sfc_partitions_module/check_partitions
  ! PURPOSE
  !   Check if the partitions are setup correctly
  ! USAGE
  !   call check_partitions()
  ! INPUTS
  !   reads module variables
  ! OUTPUT
  !   prints out info or terminates program if check fails
  ! AUTHOR
  !   L.Tong
  ! CREATION DATE
  !   2013/04/03
  ! MODIFICATION HISTORY
  ! SOURCE
  !
  subroutine check_partitions()
    ! check if the partitions are setup correctly
    use global_module,  only: ni_in_cell, iprint_init, io_lun, &
                              species_glob, numprocs, area_init
    use GenComms,       only: cq_abort, inode, ionode
    use species_module, only: nsf_species
    use memory_module,  only: reg_alloc_mem, reg_dealloc_mem, type_int

    implicit none

    integer :: tot_n_atoms, tot_n_sf, tot_n_parts, n_occupied_parts
    integer :: nsf_in_cell, n_atoms_in_parts
    integer :: ii, jj, counter, stat
    logical :: found_repeat
    integer, dimension(:), allocatable :: atoms_in_parts

    allocate(atoms_in_parts(ni_in_cell), STAT=stat)
    if (stat /= 0) then
       call cq_abort("check_partitions: Error allocating atoms_in_parts", &
                     ni_in_cell)
    end if
    call reg_alloc_mem(area_init, ni_in_cell, type_int)

    ! check total partitions
    if ((inode == ionode) .and. (iprint_init >= 2)) then
       write (io_lun, "(10x,a)", advance="no") &
            "checking if total number of partitions >= number of procs... "
    end if
    tot_n_parts = product(n_parts)
    if (tot_n_parts < numprocs) then
       call cq_abort("sfc_partitions_module: &
                      &too few partitions created for &
                      &number of processors", &
                      tot_n_parts, numprocs)
    else
       if ((inode == ionode) .and. (iprint_init >= 2)) then
          write (io_lun, "(a)") "PASS"
       end if
    end if
    ! check total occupied partitions
    if ((inode == ionode) .and. (iprint_init >= 2)) then
       write (io_lun, "(10x,a)", advance="no") &
            "checking if number of occupied partitions >= number of procs... "
    end if
    n_occupied_parts = 0
    tot_n_atoms = 0
    tot_n_sf = 0
    do ii = 1, tot_n_parts
       if (parts(ii)%n_atoms > 0) then
          n_occupied_parts = n_occupied_parts + 1
          tot_n_atoms = tot_n_atoms + parts(ii)%n_atoms
          tot_n_sf = tot_n_sf + parts(ii)%n_sf
       end if
    end do
    if (n_occupied_parts < numprocs) then
       call cq_abort("sfc_partitions_module: &
                      &too few occupied partitions for &
                      &number of processors", &
                      n_occupied_parts, numprocs)
    else
       if ((inode == ionode) .and. (iprint_init >= 2)) then
          write (io_lun, "(a)") "PASS"
       end if
    end if
    ! check total number of atoms
    if ((inode == ionode) .and. (iprint_init >= 2)) then
       write (io_lun, "(10x,a)", advance="no") &
            "checking if total N atoms in partitions == N atoms in cell... "
    end if
    if (tot_n_atoms /= ni_in_cell) then
       call cq_abort("sfc_partitions_module: &
                      &number of atoms in partitions not &
                      &equal to number of atoms in cell", &
                      tot_n_atoms, ni_in_cell)
    else
       if ((inode == ionode) .and. (iprint_init >= 2)) then
          write (io_lun, "(a)") "PASS"
       end if
    end if
    ! check total number of support functions
    if ((inode == ionode) .and. (iprint_init >= 2)) then
       write (io_lun, "(10x,a)", advance="no") &
            "checking if total N support fns in partitions == N sf in cell... "
    end if
    nsf_in_cell = 0
    do ii = 1, ni_in_cell
       nsf_in_cell = nsf_in_cell + nsf_species(species_glob(ii))
    end do
    if (tot_n_sf /= nsf_in_cell) then
       call cq_abort("sfc_partitions_module: &
                      &number of support fns in partitions not &
                      &equal to number of support fns in cell", &
                      tot_n_sf, nsf_in_cell)
    else
       if ((inode == ionode) .and. (iprint_init >= 2)) then
          write (io_lun, "(a)") "PASS"
       end if
    end if
    ! check if the atoms are assigned to partitions in a correct way
    if ((inode == ionode) .and. (iprint_init >= 2)) then
       write (io_lun, "(10x,a)", advance="no") &
            "checking if atoms are correctly assigned to partitions... "
    end if
    n_atoms_in_parts = 0
    counter = 1
    do ii = 1, tot_n_parts
       do jj = 1, parts(ii)%n_atoms
          if ((iglob_atom_part(jj,ii) > 0) .and. &
              (iglob_atom_part(jj,ii) <= ni_in_cell)) then
             n_atoms_in_parts = n_atoms_in_parts + 1
             atoms_in_parts(counter) = iglob_atom_part(jj,ii)
             counter = counter + 1
          end if
       end do
    end do
    if (n_atoms_in_parts /= ni_in_cell) then
       call cq_abort("sfc_partitions_module: &
                      &atoms are not correctly assigned to partitions.&
                      &The total number of atoms in partitions does &
                      &not equal to total number of atoms in cell")
    else
       ! check if there are duplicates in atomic indices
       found_repeat = .false.
       loop_ii: do ii = ni_in_cell, 2, -1
          do jj = 1, ii-1
             if (atoms_in_parts(jj) == atoms_in_parts(ii)) then
                found_repeat = .true.
                exit loop_ii
             end if
          end do
       end do loop_ii
       if (found_repeat) then
          deallocate(atoms_in_parts, STAT=stat)
          call cq_abort("sfc_partitions_module: &
                         &atoms are not correctly assigned to partitions.&
                         &Found duplicated atomic indices.")
       else
          if ((inode == ionode) .and. (iprint_init >= 2)) then
             write (io_lun, "(a)") "PASS"
          end if
       end if
    end if
    deallocate(atoms_in_parts, STAT=stat)
    if (stat /= 0) then
       call cq_abort("check_partitions: Error deallocating atoms_in_parts")
    end if
    call reg_dealloc_mem(area_init, ni_in_cell, type_int)

    return
  end subroutine check_partitions
  !*****


  !****f* sfc_partitions_module/print_information
  ! PURPOSE
  !   print partitioning information and statistics
  ! USAGE
  !   call print_information()
  ! INPUTS
  !   module variables
  ! OUTPUT
  !   prints information to io_lun
  ! AUTHOR
  !   L.Tong
  ! CREATION DATE
  !   2013/04/03
  ! MODIFICATION HISTORY
  ! SOURCE
  !
  subroutine print_information()

    use datatypes
    use global_module, only: iprint_init, io_lun, numprocs, area_init
    use GenComms,      only: inode, ionode, cq_abort
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_int

    implicit none

    ! local variables
    integer :: n_occupied_parts, ii, jj, ipart, stat
    integer :: max_natoms_proc,    min_natoms_proc
    integer :: max_nparts_proc,    min_nparts_proc
    integer :: max_nsf_proc,       min_nsf_proc
    integer :: max_noccparts_proc, min_noccparts_proc
    integer, dimension(:), allocatable :: n_occupied_parts_proc

    real(double) :: mean_natoms_proc, std_natoms_proc
    real(double) :: mean_nparts_proc, std_nparts_proc
    real(double) :: mean_nsf_proc, std_nsf_proc
    real(double) :: mean_noccparts_proc, std_noccparts_proc

    ! only the ionode will carry out this subroutine
    if (inode == ionode) then
       ! print system information
       select case (FSC%system_type)
       case (3)
          write (io_lun, "(/,10x,a,/)") "Detected system type: bulk"
       case (2)
          write (io_lun, "(/,10x,a,/)") "Detected system type: slab"
       case (1)
          write (io_lun, "(/,10x,a,/)") "Detected system type: chain"
       case (0)
          write (io_lun, "(/,10x,a,/)") "Detected system type: molecule"
       end select
       if (iprint_init >= 2) then
          write (io_lun, "(10x,a,3f10.6)")  &
               "Cell dimensions    (a0): ", &
               FSC%dims(1), FSC%dims(2), FSC%dims(3)
          write (io_lun, "(10x,a,2f10.6)")  &
               "System extent in x (a0): ", &
               FSC%r_atoms_min(1), FSC%r_atoms_max(1)
          write (io_lun, "(10x,a,2f10.6)")  &
               "System extent in y (a0): ", &
               FSC%r_atoms_min(2), FSC%r_atoms_max(2)
          write (io_lun, "(10x,a,2f10.6)")  &
               "System extent in z (a0): ", &
               FSC%r_atoms_min(3), FSC%r_atoms_max(3)
          write (io_lun, "(10x,a,f10.6)")   &
               "Largest gap in x   (a0): ", &
               FSC%gap(2,1) - FSC%gap(1,1)
          write (io_lun, "(10x,a,f10.6)")   &
               "Largest gap in y   (a0): ", &
               FSC%gap(2,2) - FSC%gap(1,2)
          write (io_lun, "(10x,a,f10.6)")   &
               "Largest gap in z   (a0): ", &
               FSC%gap(2,3) - FSC%gap(1,3)
          write (io_lun, "(10x,a,3l3)")     &
               "Partitioning method treating system is folded in (x,y,z): ", &
               FSC%is_folded(1), FSC%is_folded(2), FSC%is_folded(3)
       end if
       ! print partition information
       if (iprint_init >= 2) then
          write (io_lun, "(10x,a,i6)")  &
               "User requested N partitions in x: ", n_parts_user(1)
          write (io_lun, "(10x,a,i6)")  &
               "User requested N partitions in y: ", n_parts_user(2)
          write (io_lun, "(10x,a,i6)")  &
               "User requested N partitions in z: ", n_parts_user(3)
          write (io_lun, "(10x,a)")     &
               "(zero means set partitions automatically)"
          write (io_lun, "(10x,a,i6)")  &
               "Actual N partitions in x: ", n_parts(1)
          write (io_lun, "(10x,a,i6)")  &
               "Actual N partitions in y: ", n_parts(2)
          write (io_lun, "(10x,a,i6)")  &
               "Actual N partitions in z: ", n_parts(3)
          write (io_lun, "(10x,a,3i6)") &
               "Number of recursions used for 3D Hilbert curve (x,y,z): ", &
               n_divs(1), n_divs(2), n_divs(3)
          write (io_lun, "(10x,a,3f10.6)") &
               "Partition cell size (a0): ", r_part(1), r_part(2), r_part(3)
       end if
       ! print partition statistics
       if (iprint_init >= 2) then
          allocate(n_occupied_parts_proc(numprocs), STAT=stat)
          if (stat /= 0) then
             call cq_abort("print_information: &
                            &Error allocating n_occupied_parts_proc", &
                            numprocs)
          end if
          call reg_alloc_mem(area_init, numprocs, type_int)
          write (io_lun, "(10x,a)") "Partitioning Statistics:"
          write (io_lun, "(12x,a,i6)") "Total N procs: ", numprocs
          n_occupied_parts = 0
          max_nparts_proc = 0
          max_natoms_proc = 0
          max_nsf_proc = 0
          max_noccparts_proc = 0
          min_nparts_proc = huge(1)
          min_natoms_proc = huge(1)
          min_nsf_proc = huge(1)
          min_noccparts_proc = huge(1)
          mean_nparts_proc = 0.0_double
          mean_natoms_proc = 0.0_double
          mean_nsf_proc = 0.0_double
          mean_noccparts_proc = 0.0_double
          do ii = 1, numprocs
             n_occupied_parts_proc(ii) = 0
             do jj = 1, procs(ii)%n_parts
                ipart = procs(ii)%i_part_start + jj - 1
                if (parts(ipart)%n_atoms > 0) then
                   n_occupied_parts_proc(ii) = n_occupied_parts_proc(ii) + 1
                   n_occupied_parts = n_occupied_parts + 1
                end if
             end do
             mean_noccparts_proc = mean_noccparts_proc + &
                                   real(n_occupied_parts_proc(ii),double)
             mean_nparts_proc = mean_nparts_proc + &
                                real(procs(ii)%n_parts,double)
             mean_natoms_proc = mean_natoms_proc + &
                                real(procs(ii)%n_atoms,double)
             mean_nsf_proc = mean_nsf_proc + real(procs(ii)%n_sf,double)
             max_nparts_proc = max(procs(ii)%n_parts, max_nparts_proc)
             max_natoms_proc = max(procs(ii)%n_atoms, max_natoms_proc)
             max_nsf_proc = max(procs(ii)%n_sf, max_nsf_proc)
             max_noccparts_proc = max(n_occupied_parts_proc(ii), &
                                      max_noccparts_proc)
             min_nparts_proc = min(procs(ii)%n_parts, min_nparts_proc)
             min_natoms_proc = min(procs(ii)%n_atoms, min_natoms_proc)
             min_nsf_proc = min(procs(ii)%n_sf, min_nsf_proc)
             min_noccparts_proc = min(n_occupied_parts_proc(ii), &
                                      min_noccparts_proc)
          end do
          mean_nparts_proc = mean_nparts_proc / real(numprocs,double)
          mean_natoms_proc = mean_natoms_proc / real(numprocs,double)
          mean_nsf_proc = mean_nsf_proc / real(numprocs,double)
          mean_noccparts_proc = mean_noccparts_proc / real(numprocs,double)
          std_nparts_proc = 0.0_double
          std_natoms_proc = 0.0_double
          std_nsf_proc = 0.0_double
          std_noccparts_proc = 0.0_double
          do ii = 1, numprocs
             std_nparts_proc = std_nparts_proc + &
                  (real(procs(ii)%n_parts,double) - mean_nparts_proc)**2
             std_natoms_proc = std_natoms_proc + &
                  (real(procs(ii)%n_atoms,double) - mean_natoms_proc)**2
             std_nsf_proc = std_nsf_proc + &
                  (real(procs(ii)%n_sf,double) - mean_nsf_proc)**2
             std_noccparts_proc = std_noccparts_proc + &
                  (real(n_occupied_parts_proc(ii),double) - &
                   mean_noccparts_proc)**2
          end do
          std_nparts_proc = sqrt(std_nparts_proc / real(numprocs,double))
          std_natoms_proc = sqrt(std_natoms_proc / real(numprocs,double))
          std_nsf_proc = sqrt(std_nsf_proc / real(numprocs,double))
          std_noccparts_proc = sqrt(std_noccparts_proc / real(numprocs,double))
          write (io_lun, "(12x,a,i6)")      &
               "Total N partitions: ", product(n_parts)
          write (io_lun, "(12x,a,i6)")      &
               "Max   N partitions in each proc: ", max_nparts_proc
          write (io_lun, "(12x,a,i6)")      &
               "Min   N partitions in each proc: ", min_nparts_proc
          write (io_lun, "(12x,a,f10.6)")   &
               "Mean  N partitions in each proc: ", mean_nparts_proc
          write (io_lun, "(12x,a,f10.6)")   &
               "Std   N partitions in each proc: ", std_nparts_proc
          write (io_lun, "(12x,a,i6)")      &
               "Total occupied partitions: ", n_occupied_parts
          write (io_lun, "(12x,a,i6)")      &
               "Max   N occupied partitions in each proc: ", max_noccparts_proc
          write (io_lun, "(12x,a,i6)")      &
               "Min   N occupied partitions in each proc: ", min_noccparts_proc
          write (io_lun, "(12x,a,f10.6)")   &
               "Mean  N occupied partitions in each proc: ", mean_noccparts_proc
          write (io_lun, "(12x,a,f10.6)")   &
               "Std   N occupied partitions in each proc: ", std_noccparts_proc
          write (io_lun, "(12x,a,i6)")      &
               "Max   N atoms in each proc: ", max_natoms_proc
          write (io_lun, "(12x,a,i6)")      &
               "Min   N atoms in each proc: ", min_natoms_proc
          write (io_lun, "(12x,a,f10.6)")   &
               "Mean  N atoms in each proc: ", mean_natoms_proc
          write (io_lun, "(12x,a,f10.6)")   &
               "Std   N atoms in each proc: ", std_natoms_proc
          write (io_lun, "(12x,a,i6)")      &
               "Max   N support fns in each proc: ", max_nsf_proc
          write (io_lun, "(12x,a,i6)")      &
               "Min   N support fns in each proc: ", min_nsf_proc
          write (io_lun, "(12x,a,f10.6)")   &
               "Mean  N support fns in each proc: ", mean_nsf_proc
          write (io_lun, "(12x,a,f10.6,/)") &
               "Std   N support fns in each proc: ", std_nsf_proc
          deallocate(n_occupied_parts_proc, STAT=stat)
          if (stat /= 0) then
             call cq_abort("print_information: &
                            &Error deallocating n_occupied_parts_proc")
          end if
          call reg_dealloc_mem(area_init, numprocs, type_int)
       end if
    end if ! (inode == ionode)

    return
  end subroutine print_information
  !*****

end module sfc_partitions_module
