! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module mlff
! ------------------------------------------------------------------------------
! (Code area 9: general) not sure ***
! ------------------------------------------------------------------------------

!!****h* Conquest/mlff
!!  NAME
!!   mlff
!!  PURPOSE
!!   Obtain forces and energies from machine learning force field
!!  AUTHOR
!!   Jianbo Lin
!!  CREATION DATE
!!   2022/03/14
!!  MODIFICATION HISTORY
!
!!  SOURCE
!!
module mlff_type

  use datatypes
  use energy, only: disp_energy
  use GenComms,               only: cq_abort
  use timer_module,           only: start_timer,stop_timer
  use timer_stdclocks_module, only: tmr_std_allocation
  use energy, only: ml_energy

  implicit none

  save

  logical :: flag_debug_mlff = .FALSE.  ! .TRUE. ! control output level from mlff for debug
  logical :: flag_time_mlff =  .TRUE. ! .FALSE. ! control output level from mlff for timing

!! Modified from * matrix_module/matrix *
!!****s* ML_pair *
!!  NAME
!!   atom_pair
!!  PURPOSE
!!   Combines neighbour information needed to index matrices for
!!  machine learning. Save in one dimensional list.
!!  AUTHOR
!!   Jianbo Lin
!!  SOURCE
!!
  type ML_pair
    integer :: mx_nab,mx_abs
    integer :: length
    real(double)   :: rcut
    integer(integ) :: n_atoms, part_nabs ! Lengths of n_nab, i_acc

    ! index shift of neighbours
    ! No. of neighbours of each primary atom i
    integer(integ),pointer :: n_nab(:)
    ! accumulator for no. of neighbours of primary atom i
    integer(integ),pointer :: i_acc(:)
    integer(integ),pointer :: i_species(:)  ! species of target atom i

    ! useful information for machine learning
    ! the information from each neighbour
    integer(integ),pointer :: j_species(:)  ! species of neighbour atom j
    ! save vector of rij and r, or calcuate directly later
    real(double),pointer :: rijx(:)
    real(double),pointer :: rijy(:)
    real(double),pointer :: rijz(:)
    real(double),pointer :: radius(:)
  end type ML_pair
!!***
  type matrix_ML
    integer :: mx_nab,mx_abs
    integer :: length
    ! Useful information for machine learning
    real(double)   :: rcut
    integer(integ) :: array_posn ! Starting position in array
    integer(integ) :: offset     ! This plus i_acc gives us i_dir
    integer(integ) :: nd_offset  ! This plus i_acc gives us i_dir
    integer(integ) :: n_atoms, part_nabs, part_nd_nabs ! Lengths of n_nab, i_acc, i_nd_acc
    integer(integ),pointer :: onsite(:)
    ! From here on, we point to the comms array
    integer(integ),pointer :: n_nab(:)     ! Neighbours of primary atom i
    integer(integ),pointer :: i_acc(:)     ! accumulator for atom in neigh list  : ist_j = amat(nn)%i_acc(i)+j-1
    integer(integ),pointer :: i_nd_acc(:)  ! accumulator for SFs in neigh list
    integer(integ),pointer :: npxyz(:)     ! offset of CS part of atom in n. list

    ! Information from each neighbour
    integer(integ),pointer :: i_species(:)  ! species of target atom i, dim: natom
    ! Save vector of rij and r, or calcuate directly later
    integer(integ),pointer :: j_species(:)  ! species of neighbour atom j  {[ith naba], [i+1th naba],}, dim
    integer(integ),pointer :: i_part(:)     ! partition no. of atom in neigh list
    integer(integ),pointer :: i_seq(:)      ! part-seq no. of atom in neigh list
    real(double),pointer :: dx(:)
    real(double),pointer :: dy(:)
    real(double),pointer :: dz(:)
    real(double),pointer :: radius(:)

    ! Information of three body
    integer(integ),pointer :: triplet_acc(:) ! accumulator for atom in triplet list
    integer(integ),pointer :: n_triplet(:)   ! Three body triplets of primary atom i
    integer(integ),pointer :: ist_triplet(:) ! seq no. of atom in neigh (pair) list
    integer(integ) :: part_triplets           ! no. of triplets
  end type matrix_ML

  !! Allow each atom to have different dimensions of feature
  type features_atomic
    integer(integ) :: dim_feature ! dimension_feature
    real(double), allocatable  :: fpx(:) ! feature for atomic force in x direction of natoms
    real(double), allocatable  :: fpy(:) ! feature for atomic force in x direction of natoms
    real(double), allocatable  :: fpz(:) ! feature for atomic force in x direction of natoms
    real(double), allocatable  :: fp(:)  ! feature for atomic energy of natoms
  end type features_atomic

  type features_ML
    integer(integ) :: n_atoms ! Lengths of n_nab,
    type(features_atomic),  allocatable :: id_atom(:) ! allow each atom to have different dimensions of feature
  end type features_ML

  type species_order
    integer, dimension(:), allocatable :: d2
    integer, dimension(:,:), allocatable :: d3
  end type species_order

!!****s* ML_type/ML_acsf *
!!  NAME
!!   ML_acsf - machine learning descriptor: atomic centered symmetry function
!!  PURPOSE
!!   Defines the parameters for acsf functions
!!  AUTHOR
!!   Jianbo Lin
!!  CREATION DATE
!!   2022/03/14
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  type g2_param
    real(double), allocatable :: eta(:)
    real(double), allocatable :: rs(:)
    real(double), allocatable :: scales(:)
    real(double), allocatable :: coefs(:)
  end type g2_param

  type g3_param
    real(double), pointer :: kappa(:)
    real(double), pointer :: scales(:)
    real(double), pointer :: coefs(:)
  end type g3_param

  type g45_param
    real(double), pointer :: eta(:)
    real(double), pointer :: zeta(:)
    real(double), pointer :: lambda(:)
    real(double), pointer :: scales(:)
    real(double), pointer :: coefs(:)
  end type g45_param

  type acsf_param
    ! Number of g2,g3,g4,g5
    integer :: num_g2,num_g3,num_g4,num_g5                    ! total number of two-body or three-body terms
    integer, dimension(:), allocatable :: nums_g2,nums_g3,nums_g4,nums_g5   ! numbers of two-body and three-body terms
    integer, dimension(:), allocatable :: nums_g2_acc,nums_g3_acc,nums_g4_acc,nums_g5_acc  ! start index of each component

    integer :: n_species                             ! number of species in this descriptor
    character(len=2), allocatable :: species_lst(:)  ! dimension n_species
    type(species_order)  :: species_orders           ! order of two-body and three-body terms in descriptor
    ! Cutoff
    real(double) :: rcut
    ! Arrays
    type(g2_param),  allocatable  :: params_g2(:)    ! two lists of Eta and Rs parameters for G2 functions
    type(g3_param),  allocatable  :: params_g3(:)    ! a list of kapa parameters for G3 functions
    type(g45_param),  allocatable :: params_g4(:)    ! three lists of Eta, zeta and lamda parameters for G4 functions
    type(g45_param),  allocatable :: params_g5(:)    ! three lists of Eta, zeta and lamda parameters for G5 functions
    !Model parameter in one dimension
    integer :: dim_coef
    real(double),allocatable :: coef(:) ! this coef here is standardscale * coefficent of Linear-Regression
  end type acsf_param
!!***

!!****s* ML_type/ML_split *
!!  NAME
!!   ML_split - machine learning descriptor: split type symmetry function
!!  PURPOSE
!!   Defines the parameters for split functions
!!  AUTHOR
!!   Jianbo Lin
!!  CREATION DATE
!!   2022/03/14
!!  MODIFICATION HISTORY
!!  SOURCE
!!
  type b3_param
    real(double), allocatable :: eta1(:)
    real(double), allocatable :: eta2(:)
    real(double), allocatable :: eta3(:)
    real(double), allocatable :: scales(:)
    real(double), allocatable :: coefs(:)
  end type b3_param

  type split_param
    ! Number of 2b,3b
    integer          :: num_2b,num_3b                  ! total number of two-body or three-body terms
    integer, dimension(:), allocatable :: nums_2b,nums_3b            ! numbers of two-body and three-body terms
    integer, dimension(:), allocatable :: nums_2b_acc, nums_3b_acc    ! start index of each component

    integer :: n_species                             ! number of species in this descriptor
    character(len=2), allocatable :: species_lst(:)  ! dimension n_species
    type(species_order)  :: species_orders           ! order of two-body and three-body terms in descriptor
    ! Cutoff
    real(double) :: rcut
    ! Arrays
    type(g2_param), allocatable :: params_2b(:)    ! a list of Eta and Rs parameters for two-body functions
    type(b3_param), allocatable :: params_3b(:)    ! a list of Eta1, Eta2 and Eta3 parameters for three-body functions
    !Model parameter in one dimension
    integer :: dim_coef                 ! dimension of feature
    real(double),allocatable :: coef(:) ! this coef here is standardscale * coefficent of Linear-Regression
  end type split_param

  type descriptor_param
    character(len=80)          :: descriptor_type
    type(acsf_param),  allocatable :: acsf(:)
    type(split_param), allocatable :: split(:)
    integer,allocatable :: dim_coef_lst(:) ! collect dimension of feature for each species
  end type descriptor_param

  type(descriptor_param) :: params_ML
  !integer, dimension(:,:), allocatable :: species_order
!!***

  interface swap ! swap two data
     module procedure int_swap
     module procedure double_swap
  end interface

  interface read_ml ! read machine learning model
     module procedure read_acsf2b
     module procedure read_split
  end interface

contains


!!****s* ML_type/set_ML *
!!  NAME
!!   set_ML
!!  PURPOSE
!!   Read the parameters from machine learning infomation file
!!  AUTHOR
!!   Jianbo Lin
!!  CREATION DATE
!!   2022/03/14
!!  MODIFICATION HISTORY
!!   2023/08/22 J.Lin
!!   Added checking if file exists
!!  SOURCE
!!
  subroutine set_ML

    use GenComms, ONLY: inode, ionode

    implicit none
    integer                    :: file_id=2023, stat, n_species,ii
    character(len=80)          :: filename
    character(len=80)          :: descriptor_type, comment
    ! Todo: filename
    filename='ml_pot.txt'
    ! After we have neighbor information and descriptor information
    ! read and set
    open(unit=file_id,file=filename,iostat=stat,status='old')
    if(stat/=0) &
      call cq_abort('open ml_pot.txt : error no such file')
    read(file_id,*) comment
    read(file_id,*) descriptor_type
    close(file_id)
    descriptor_type=TRIM(descriptor_type)
    params_ML%descriptor_type = descriptor_type

    !! TODO: check flag_descriptortype
    if (inode == ionode) then
      write(*,*) 'flag_descriptortype:', descriptor_type
    end if

    if (descriptor_type == 'split2b3b' .or. descriptor_type == 'Split2b3b') then
      call read_ml(filename, params_ML%split)
      n_species=params_ML%split(1)%n_species
      allocate(params_ML%dim_coef_lst(n_species),STAT=stat)
      if(stat/=0) &
        call cq_abort('params_ML%dim_coef_lst : error allocating memory to n_species')
      do ii=1, n_species
        params_ML%dim_coef_lst(ii) = params_ML%split(ii)%dim_coef
      end do ! ii1
    elseif (descriptor_type == 'acsf2b' .or. descriptor_type == 'Acsf2b') then
      call read_ml(filename, params_ML%acsf)
      n_species=params_ML%acsf(1)%n_species
      allocate(params_ML%dim_coef_lst(n_species),STAT=stat)
      if(stat/=0) &
        call cq_abort('params_ML%dim_coef_lst : error allocating memory to n_species')
      do ii=1, n_species
        params_ML%dim_coef_lst(ii) = params_ML%acsf(ii)%dim_coef
      end do ! ii1
    end if

    if (inode == ionode) then
      write(*,*) 'descriptor dimensions:'
      write(*,*) params_ML%dim_coef_lst
    end if
  end subroutine set_ML

!!****f* matrix_module/allocate_matrix_ML *
!!
!!  NAME
!!   allocate_matrix_ML
!!  USAGE
!!
!!  PURPOSE
!!   Allocate memory to matrix derived type
!!  INPUTS
!!
!!
!!  USES
!!
!!  AUTHOR
!!   Jianbo Lin
!!  CREATION DATE
!!   2023/07/20
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine allocate_matrix_ML(mat_ML,part_on_node,mx_part)
    use matrix_module, ONLY: matrix
    use GenComms, ONLY: cq_abort
    implicit none

    ! Passed variables
    integer :: mx_part,part_on_node
    type(matrix_ML), dimension(:)  :: mat_ML

    ! Local variables
    integer :: mx_naba, i, stat

    call start_timer(tmr_std_allocation)
    do i=1,part_on_node
      ! Dimensions of neighbors
      mx_naba = mx_part*mat_ML(i)%mx_abs
      if (mx_naba .gt. 0) then
        allocate(mat_ML(i)%i_seq(mx_naba),STAT=stat)
        if(stat/=0) &
          call cq_abort('alloc_mat_ML: error allocating memory to i_seq')
        allocate(mat_ML(i)%i_part(mx_naba),STAT=stat)
        if(stat/=0) &
          call cq_abort('alloc_mat_ML: error allocating memory to i_part')
        allocate(mat_ML(i)%radius(mx_naba),STAT=stat)
        if(stat/=0) &
          call cq_abort('alloc_mat_ML: error allocating memory to radius')
        allocate(mat_ML(i)%j_species(mx_naba),STAT=stat)
        if(stat/=0) &
          call cq_abort('alloc_mat_ML: error allocating memory to j_species')
        allocate(mat_ML(i)%dx(mx_naba),STAT=stat)
        if(stat/=0) &
          call cq_abort('alloc_mat_ML: error allocating memory to dx')
        allocate(mat_ML(i)%dy(mx_naba),STAT=stat)
        if(stat/=0) &
          call cq_abort('alloc_mat_ML: error allocating memory to dy')
        allocate(mat_ML(i)%dz(mx_naba),STAT=stat)
        if(stat/=0) &
          call cq_abort('alloc_mat_ML: error allocating memory to dz')
      else
        write(*, *) 'Warning: allocate matrix ML, partition, mx_naba', i, mx_naba
      end if ! mx_naba
      ! Dimensions of atoms in processes
      if (mx_part .gt. 0) then
        allocate(mat_ML(i)%onsite(mx_part),STAT=stat)
        if(stat/=0) &
          call cq_abort('alloc_mat_ML: error allocating memory to onsite')
        allocate(mat_ML(i)%n_nab(mx_part),STAT=stat)
        if(stat/=0) &
          call cq_abort('alloc_mat_ML: error allocating memory to n_nab')
        allocate(mat_ML(i)%i_acc(mx_part),STAT=stat)
        if(stat/=0) &
          call cq_abort('alloc_mat_ML: error allocating memory to i_acc')
        allocate(mat_ML(i)%i_nd_acc(mx_part),STAT=stat)
        if(stat/=0) &
          call cq_abort('alloc_mat_ML: error allocating memory to i_nd_acc')
        allocate(mat_ML(i)%i_species(mx_part),STAT=stat)
        if(stat/=0) &
          call cq_abort('alloc_mat_ML: error allocating memory to i_species')
      else
        write(*, *) 'Warning: allocate matrix ML, partition, mx_part', i, mx_part
      end if ! mx_part
    enddo ! i, part_on_node
    call stop_timer(tmr_std_allocation)
    return
  end subroutine allocate_matrix_ML
!!***

!!****f* matrix_module/deallocate_matrix_ML *
!!
!!  NAME
!!   deallocate_matrix_ML
!!  USAGE
!!
!!  PURPOSE
!!   Deallocate memory to matrix derived type
!!  INPUTS
!!
!!
!!  USES
!!
!!  AUTHOR
!!   Jianbo Lin
!!  CREATION DATE
!!   2023/07/20
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine deallocate_matrix_ML(mat_ML,part_on_node,mx_part)
    use matrix_module, ONLY: matrix
    use GenComms, ONLY: cq_abort
    implicit none

    ! Passed variables
    integer :: mx_part,part_on_node
    type(matrix_ML), dimension(:)  :: mat_ML

    ! Local variables
    integer :: mx_naba, i, stat

    if(index(params_ML%descriptor_type, '3b') .gt. 0)then
      call deallocate_matrix_ML_triplet(mat_ML, part_on_node)
    end if
    call start_timer(tmr_std_allocation)
    do i=part_on_node,1,-1
      ! Dimensions of neighbors
      mx_naba = mx_part*mat_ML(i)%mx_abs
      if (mx_naba .gt. 0) then
        deallocate(mat_ML(i)%i_seq,STAT=stat)
        if(stat/=0) &
          call cq_abort('dealloc_mat_ML: error deallocating memory to i_seq')
        deallocate(mat_ML(i)%i_part,STAT=stat)
        if(stat/=0) &
          call cq_abort('dealloc_mat_ML: error deallocating memory to i_part')
        deallocate(mat_ML(i)%radius,STAT=stat)
        if(stat/=0) &
          call cq_abort('dealloc_mat_ML: error deallocating memory to radius')
        deallocate(mat_ML(i)%j_species,STAT=stat)
        if(stat/=0) &
          call cq_abort('dealloc_mat_ML: error deallocating memory to j_species')
        deallocate(mat_ML(i)%dx,STAT=stat)
        if(stat/=0) &
          call cq_abort('dealloc_mat_ML: error deallocating memory to dx')
        deallocate(mat_ML(i)%dy,STAT=stat)
        if(stat/=0) &
          call cq_abort('dealloc_mat_ML: error deallocating memory to dy')
        deallocate(mat_ML(i)%dz,STAT=stat)
        if(stat/=0) &
          call cq_abort('dealloc_mat_ML: error deallocating memory to dz')
      else
          write(*, *) 'Warning: deallocate matrix ML, partition, mx_naba', i, mx_naba
      end if ! mx_naba
      ! Dimensions of atoms in processes
      if (mx_part .gt. 0) then
         deallocate(mat_ML(i)%onsite,STAT=stat)
         if(stat/=0) &
            call cq_abort('dealloc_mat_ML: error deallocating memory to onsite')
         deallocate(mat_ML(i)%n_nab,STAT=stat)
         if(stat/=0) &
            call cq_abort('dealloc_mat_ML: error deallocating memory to n_nab')
         deallocate(mat_ML(i)%i_acc,STAT=stat)
         if(stat/=0) &
            call cq_abort('dealloc_mat_ML: error deallocating memory to i_acc')
         deallocate(mat_ML(i)%i_nd_acc,STAT=stat)
         if(stat/=0) &
            call cq_abort('dealloc_mat_ML: error deallocating memory to i_nd_acc')
         deallocate(mat_ML(i)%i_species,STAT=stat)
         if(stat/=0) &
            call cq_abort('dealloc_mat_ML: error deallocating memory to i_species')
      else
          write(*, *) 'Warning: deallocate matrix ML, partition, mx_part', i, mx_part
      end if ! mx_part
    enddo ! i, part_on_node
    call stop_timer(tmr_std_allocation)
    return
  end subroutine deallocate_matrix_ML
!!***


!!****f* matrix_module/allocate_matrix_ML_triplet *
!!
!!  NAME
!!   allocate_matrix_ML_triplet
!!  USAGE
!!
!!  PURPOSE
!!   Allocate memory to matrix derived type
!!  INPUTS
!!
!!
!!  USES
!!
!!  AUTHOR
!!   Jianbo Lin
!!  CREATION DATE
!!   2023/07/22
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine allocate_matrix_ML_triplet(mat_ML,part_on_node)
    use matrix_module, ONLY: matrix
    use GenComms, ONLY: cq_abort
    implicit none

    ! Passed variables
    integer :: part_on_node
    type(matrix_ML), dimension(:)  :: mat_ML

    ! Local variables
    integer :: part_nabs, part_nabs_sum,nn,ni, stat

    call start_timer(tmr_std_allocation)
    do nn=1,part_on_node
      if (mat_ML(nn)%n_atoms .gt. 0) then
        ! Total dimensions of neighbors
        part_nabs = mat_ML(nn)%part_nabs
        allocate(mat_ML(nn)%triplet_acc(part_nabs),STAT=stat)
        if(stat/=0) &
          call cq_abort('alloc_mat_ML_triplet: error allocating memory to triplet_acc')
        allocate(mat_ML(nn)%n_triplet(part_nabs),STAT=stat)
        if(stat/=0) &
          call cq_abort('alloc_mat_ML_triplet: error allocating memory to n_triplet')

        part_nabs_sum = 0
        do ni=1, mat_ML(nn)%n_atoms
          part_nabs_sum = part_nabs_sum + mat_ML(nn)%n_nab(ni) **2
        end do
        allocate(mat_ML(nn)%ist_triplet(part_nabs_sum),STAT=stat)
        if(stat/=0) &
          call cq_abort('alloc_mat_ML_triplet: error allocating memory to ist_triplet')
      else
         write(*, *) 'Warning: alloc_mat_ML_triplet, no atom in partition,', nn, mat_ML(nn)%n_atoms
      end if
    enddo ! i, part_on_node
    call stop_timer(tmr_std_allocation)
    return
  end subroutine allocate_matrix_ML_triplet
!!***


!!****f* matrix_module/deallocate_matrix_ML_triplet *
!!
!!  NAME
!!   deallocate_matrix_ML_triplet
!!  USAGE
!!
!!  PURPOSE
!!   Deallocate memory to matrix derived type
!!  INPUTS
!!
!!
!!  USES
!!
!!  AUTHOR
!!   Jianbo Lin
!!  CREATION DATE
!!   2023/07/22
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine deallocate_matrix_ML_triplet(mat_ML,part_on_node)
    use matrix_module, ONLY: matrix
    use GenComms, ONLY: cq_abort
    implicit none

    ! Passed variables
    integer :: part_on_node
    type(matrix_ML), dimension(:)  :: mat_ML

    ! Local variables
    integer :: part_nabs, part_nabs_sum,nn,stat

    call start_timer(tmr_std_allocation)
    do nn=part_on_node, 1, -1
      ! Total dimensions of neighbors
      part_nabs = mat_ML(nn)%part_nabs
      if (part_nabs .gt. 0) then
        deallocate(mat_ML(nn)%triplet_acc,STAT=stat)
        if(stat/=0) &
          call cq_abort('dealloc_mat_ML_triplet: error deallocating memory to triplet_acc')
        deallocate(mat_ML(nn)%n_triplet,STAT=stat)
        if(stat/=0) &
          call cq_abort('dealloc_mat_ML_triplet: error deallocating memory to n_triplet')
        deallocate(mat_ML(nn)%ist_triplet,STAT=stat)
        if(stat/=0) &
          call cq_abort('dealloc_mat_ML_triplet: error deallocating memory to ist_triplet')
      else
        write(*, *) 'Warning: deallocate matrix ML triplet: partition, part_nabs', nn, part_nabs
      end if
    enddo ! i, part_on_node
    call stop_timer(tmr_std_allocation)
    return
  end subroutine deallocate_matrix_ML_triplet
!!***

!!****f* matrix_module/ini_species_order *
  subroutine ini_species_order(species_orders, n_species)
    implicit none

    ! Passed variables
    type(species_order) :: species_orders
    integer                 :: n_species

    ! Local variables
    integer :: i, j, k, shift,stat

    allocate(species_orders%d2(n_species),STAT=stat)
    if(stat/=0) then
        call cq_abort('species_orders d2: error allocating memory to n_species')
    endif
    allocate(species_orders%d3(n_species, n_species),STAT=stat)
    if(stat/=0) then
        call cq_abort('species_orders d3: error allocating memory to (n_species, n_species)')
    endif
    shift=0
    do i=1, n_species
      species_orders%d2(i) = i
      do j=1, n_species
        if (j >= i) then
          species_orders%d3(i,j) = shift + j -i +1
        else
          species_orders%d3(i,j) = species_orders%d3(j,i)
        end if
      end do
      shift = shift + n_species - i +1
    end do
    return
  end subroutine ini_species_order
!!***

!!****f* matrix_module/ini_species_order *
  subroutine ini_species_order_center(species_orders, n_species, i_species)
    implicit none

    ! Passed variables
    type(species_order) :: species_orders
    integer                 :: n_species, i_species

    ! Local variables
    integer :: i, j, k, shift,stat
    integer :: i_d2, i_d3

    allocate(species_orders%d2(n_species),STAT=stat)
    if(stat/=0) then
        call cq_abort('species_orders d2: error allocating memory to n_species')
    endif
    allocate(species_orders%d3(n_species, n_species),STAT=stat)
    if(stat/=0) then
        call cq_abort('species_orders d3: error allocating memory to (n_species, n_species)')
    endif

    !! center atom to be the first order
    !! X-X, and X-XX
    species_orders%d2(i_species) = 1
    species_orders%d3(i_species,i_species) = 1
    i_d2=2
    i_d3=2
    !! X-XA,X-XB,...X-AX
    do i=1, n_species
      do j=1, n_species
        if (i == i_species .or. j==i_species) then
            if (j > i ) then
               species_orders%d3(i,j) = i_d3
               i_d3 = i_d3 + 1
            else
               species_orders%d3(i,j) = species_orders%d3(j,i)
            end if
        end if
      end do
    end do

    ! X-A, and X-AB, X-BA
    do i=1, n_species
       if (i /= i_species) then
          species_orders%d2(i) = i_d2
          i_d2 = i_d2 + 1
       end if
       do j=1, n_species
        if (i /= i_species .and. j/=i_species) then
            if (j >= i ) then
               species_orders%d3(i,j) = i_d3
               i_d3 = i_d3 + 1
            else
               species_orders%d3(i,j) = species_orders%d3(j,i)
            end if
        end if
       end do
    end do
    return
  end subroutine ini_species_order_center
!!***


!!****f* matrix_module/allocate_feature_ML *
!!
!!  NAME
!!   allocate_feature_ML
!!  USAGE
!!
!!  PURPOSE
!!   A matrix to collect information of feature for each partition.
!!  INPUTS
!!
!!
!!  USES
!!
!!  AUTHOR
!!   Jianbo Lin
!!  CREATION DATE
!!   2023/04/22
!!  MODIFICATION HISTORY
!!   2023/08/24 J.Lin
!!   Available different feature dimensions for each type of atom
!!  SOURCE
!!
  subroutine allocate_features_ML(amat_ML, amat_features_ML,prim,feature_dim_lst)
    use basic_types, ONLY: primary_set
    use GenComms, ONLY: cq_abort, inode, ionode
    implicit none

    ! Passed variables
    type(matrix_ML), dimension(:)    :: amat_ML
    type(features_ML), dimension(:)  :: amat_features_ML
    type(primary_set), intent(in)    :: prim
    integer, dimension(:)            :: feature_dim_lst

    ! Local variables
    integer ::  nn, ii, i_species, n_atoms, feature_dim, stat

    call start_timer(tmr_std_allocation)
    do nn=1,prim%groups_on_node ! Partitions in primary set
      !pair check!write(*,*) 'nn loop, inode=',inode
      if(prim%nm_nodgroup(nn).gt.0) then  ! Are there atoms ?
        if (amat_ML(nn)%n_atoms .ne. prim%nm_nodgroup(nn)) then ! Redundant, but useful
          call cq_abort('amat(nn)%n_atoms not equal prim%nm_nodgroup(nn)')
        end if
        !write(io_lun,*) 'Starting group with atoms: ',nn,prim%nm_nodgroup(nn)
        n_atoms = amat_ML(nn)%n_atoms
        allocate(amat_features_ML(nn)%id_atom(n_atoms),STAT=stat)
        if(stat/=0) &
          call cq_abort('alloc_mat_features_ML: error allocating memory to id_atom')

        do ii=1,prim%nm_nodgroup(nn)  ! Loop over atoms in partition
          i_species = amat_ML(nn)%i_species(ii)
          feature_dim = feature_dim_lst(i_species)
          allocate(amat_features_ML(nn)%id_atom(ii)%fpx(feature_dim),STAT=stat)
          if(stat/=0) &
            call cq_abort('alloc_mat_features_ML: error allocating memory to fpx')
          allocate(amat_features_ML(nn)%id_atom(ii)%fpy(feature_dim),STAT=stat)
          if(stat/=0) &
            call cq_abort('alloc_mat_features_ML: error allocating memory to fpy')
          allocate(amat_features_ML(nn)%id_atom(ii)%fpz(feature_dim),STAT=stat)
          if(stat/=0) &
            call cq_abort('alloc_mat_features_ML: error allocating memory to fpz')
          allocate(amat_features_ML(nn)%id_atom(ii)%fp(feature_dim),STAT=stat)
          if(stat/=0) &
            call cq_abort('alloc_mat_features_ML: error allocating memory to fp')

          amat_features_ML(nn)%id_atom(ii)%fpx = 0.0
          amat_features_ML(nn)%id_atom(ii)%fpy = 0.0
          amat_features_ML(nn)%id_atom(ii)%fpz = 0.0
          amat_features_ML(nn)%id_atom(ii)%fp = 0.0
          end do ! ii1
      end if
    end do ! nn1
    call stop_timer(tmr_std_allocation)
    return
  end subroutine allocate_features_ML
!!***

!!****f* matrix_module/deallocate_feature_ML *
!!
!!  NAME
!!   deallocate_feature_ML
!!  USAGE
!!
!!  PURPOSE
!!   Deallocate the matrix for collecting information of feature.
!!  INPUTS
!!
!!
!!  USES
!!
!!  AUTHOR
!!   Jianbo Lin
!!  CREATION DATE
!!   2023/04/22
!!  MODIFICATION HISTORY
!!   2023/08/24 J.Lin
!!   Available different feature dimensions for each type of atom
!!  SOURCE
!!
  subroutine deallocate_features_ML(amat_ML,amat_features_ML,part_on_node)
    use matrix_module, ONLY: matrix
    use GenComms, ONLY: cq_abort, inode, ionode
    implicit none

    ! Passed variables
    integer :: part_on_node
    type(matrix_ML), dimension(:)  :: amat_ML
    type(features_ML), dimension(:) :: amat_features_ML

    ! Local variables
    integer ::  nn,ii,n_atoms,stat

    call start_timer(tmr_std_allocation)
    do nn=part_on_node,1,-1
      ! n_atoms in this partition
      n_atoms = amat_ML(nn)%n_atoms
      if (n_atoms .gt.0) then
        do ii = n_atoms, 1, -1
          deallocate(amat_features_ML(nn)%id_atom(ii)%fp,STAT=stat)
          if(stat/=0) &
            call cq_abort('dealloc_mat_features_ML: error deallocating memory to fp')
          deallocate(amat_features_ML(nn)%id_atom(ii)%fpz,STAT=stat)
          if(stat/=0) &
            call cq_abort('dealloc_mat_features_ML: error deallocating memory to fpz')
          deallocate(amat_features_ML(nn)%id_atom(ii)%fpy,STAT=stat)
          if(stat/=0) &
            call cq_abort('dealloc_mat_features_ML: error deallocating memory to fpy')
          deallocate(amat_features_ML(nn)%id_atom(ii)%fpx,STAT=stat)
          if(stat/=0) &
            call cq_abort('dealloc_mat_features_ML: error deallocating memory to fpx')
        end do ! ii1
        deallocate(amat_features_ML(nn)%id_atom,STAT=stat)
        if(stat/=0) &
          call cq_abort('dealloc_mat_features_ML: error deallocating memory to id_atom')
      else
        if(flag_debug_mlff) &
          write(*,*) 'Warning: No atoms in this partition deallocate_features_ML', inode, n_atoms
      end if
    enddo
    call stop_timer(tmr_std_allocation)
    return
  end subroutine deallocate_features_ML
!!***


!!****f* matrix_module/allocate_b2_param *
  subroutine allocate_b2_param(params, num_dim)
    implicit none

    ! Passed variables
    type(g2_param)          :: params
    integer                 :: num_dim

    ! Local variables
    integer                 :: stat

    allocate(params%eta(num_dim),STAT=stat)
    if(stat/=0) &
      call cq_abort('alloc_b2_param: error allocating memory to eta')
    allocate(params%rs(num_dim),STAT=stat)
    if(stat/=0) &
      call cq_abort('alloc_b2_param: error allocating memory to rs')
    allocate(params%scales(num_dim),STAT=stat)
    if(stat/=0) &
      call cq_abort('alloc_b2_param: error allocating memory to scales')
    allocate(params%coefs(num_dim),STAT=stat)
    if(stat/=0) &
      call cq_abort('alloc_b2_param: error allocating memory to coefs')
    return
  end subroutine allocate_b2_param
!!***

!!****f* matrix_module/deallocate_b2_param *
  subroutine deallocate_b2_param(params, num_b2)
    implicit none

    ! Passed variables
    type(g2_param), dimension(:)  :: params
    integer                       :: num_b2

    ! Local variables
    integer                       :: i,stat

    call start_timer(tmr_std_allocation)
    do i=num_b2, 1, -1
      deallocate(params(i)%eta,STAT=stat)
      if(stat/=0) &
        call cq_abort('dealloc_b2_param: error deallocating memory to eta')
      deallocate(params(i)%rs,STAT=stat)
      if(stat/=0) &
        call cq_abort('dealloc_b2_param: error deallocating memory to rs')
      deallocate(params(i)%scales,STAT=stat)
      if(stat/=0) &
        call cq_abort('dealloc_b2_param: error deallocating memory to scales')
      deallocate(params(i)%coefs,STAT=stat)
      if(stat/=0) &
        call cq_abort('dealloc_b2_param: error deallocating memory to coefs')
    end do
    call stop_timer(tmr_std_allocation)
    return
  end subroutine deallocate_b2_param
!!***

!!****f* matrix_module/allocate_b3_param *
  subroutine allocate_b3_param(params, num_dim)
    implicit none

    ! Passed variables
    type(b3_param)          :: params
    integer                 :: num_dim

    ! Local variables
    integer                 :: stat

    allocate(params%eta1(num_dim),STAT=stat)
    if(stat/=0) &
      call cq_abort('alloc_b3_param: error allocating memory to eta1')
    allocate(params%eta2(num_dim),STAT=stat)
    if(stat/=0) &
      call cq_abort('alloc_b3_param: error allocating memory to eta2')
    allocate(params%eta3(num_dim),STAT=stat)
    if(stat/=0) &
      call cq_abort('alloc_b3_param: error allocating memory to eta3')
    allocate(params%scales(num_dim),STAT=stat)
    if(stat/=0) &
      call cq_abort('alloc_b3_param: error allocating memory to scales')
    allocate(params%coefs(num_dim),STAT=stat)
    if(stat/=0) &
      call cq_abort('alloc_b3_param: error allocating memory to coefs')
    return
  end subroutine allocate_b3_param
!!***

!!****f* matrix_module/deallocate_b3_param *
  subroutine deallocate_b3_param(params, num_b3)
    implicit none

    ! Passed variables
    type(b3_param), dimension(:), allocatable   :: params
    integer                                     :: num_b3

    ! Local variables
    integer                                     :: i,stat

    call start_timer(tmr_std_allocation)
    do i=num_b3, 1, -1
      deallocate(params(i)%eta1,STAT=stat)
      if(stat/=0) &
        call cq_abort('dealloc_b3_param: error deallocating memory to eta1')
      deallocate(params(i)%eta2,STAT=stat)
      if(stat/=0) &
        call cq_abort('dealloc_b3_param: error deallocating memory to eta2')
      deallocate(params(i)%eta3,STAT=stat)
      if(stat/=0) &
        call cq_abort('dealloc_b3_param: error deallocating memory to eta3')
      deallocate(params(i)%scales,STAT=stat)
      if(stat/=0) &
        call cq_abort('dealloc_b3_param: error deallocating memory to scales')
      deallocate(params(i)%coefs,STAT=stat)
      if(stat/=0) &
        call cq_abort('dealloc_b3_param: error deallocating memory to coefs')
    end do
    call stop_timer(tmr_std_allocation)
    return
  end subroutine deallocate_b3_param
!!***

!!****f* matrix_module/get_b2_param *
  subroutine get_b2_param(file_id, params, dims)
    implicit none

    ! Passed variables
    integer                 :: file_id
    type(g2_param)          :: params
    integer                 :: dims

    ! Local variables
    integer                 :: i


    do i =1, dims
      read(file_id,*) params%eta(i),params%rs(i),params%scales(i),params%coefs(i)
    end do
    !read(file_id,*) params%eta
    !read(file_id,*) params%rs
    !read(file_id,*) params%scales
    !read(file_id,*) params%coefs
    return
  end subroutine get_b2_param
!!***

!!****f* matrix_module/get_b3_param *
  subroutine get_b3_param(file_id, params, dims)
    implicit none

    ! Passed variables
    integer                 :: file_id
    type(b3_param)          :: params
    integer                 :: dims

    ! Local variables
    integer                 :: i

    do i =1, dims
      read(file_id,*) params%eta1(i),params%eta2(i),params%eta3(i),params%scales(i),params%coefs(i)
    end do
    return
  end subroutine get_b3_param
!!***


!!****f*  ML_type/read_acsf2b *
!!
!!  NAME
!!   read_acsf2b
!!  PURPOSE
!!
!!  AUTHOR
!!   Jianbo Lin
!!  CREATION DATE
!!   2022/03/14
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine read_acsf2b(filename, descriptor_params)

    use GenComms, ONLY: cq_abort, inode, ionode,my_barrier
    use global_module, ONLY: area_integn
    use memory_module, ONLY: reg_alloc_mem, type_int

    implicit none

    ! Passed variables
    type(acsf_param),allocatable  :: descriptor_params(:)
    character(len=*)           :: filename

    ! Local variables
    integer                    :: ii, jj, kk, dim_2b, n_species, n_2b_pair
    integer                    :: a,b,shift, num_g2
    integer                    :: stat
    integer                    :: file_id=2022
    integer, allocatable       :: dims_2b(:)
    character(len=180)         :: comment
    character(len=30)          :: eta_type, delimiter

    delimiter=' '
    open(unit=file_id,file=filename,iostat=stat,status='old')
    read(file_id,*, end=999) comment
    read(file_id,*, end=999) comment
    read(file_id,*, end=999) comment
    read(file_id,*) n_species
    !! TODO: remove in the future about eta type
    read(file_id,*, end=999) comment
    read(file_id,*) eta_type
    eta_type=TRIM(eta_type)

    !! Allocate n_spieces of descriptor_params
    allocate(descriptor_params(n_species),STAT=stat)
    if(stat/=0) &
      call cq_abort('acsf2b descriptor_params: error allocating memory to n_species')

    do ii=1, n_species
      ! species_orders only depending on input of species order
      !call ini_species_order(descriptor_params(ii)%species_orders, n_species)
      ! species_orders depending on center atom according to the input of species order
      call ini_species_order_center(descriptor_params(ii)%species_orders, n_species, ii)
      if (inode == ionode) then
        write(*,*) 'read_split2b3b: species_order d2 ', descriptor_params(ii)%species_orders%d2
        write(*,*) 'read_split2b3b: species_order d3', descriptor_params(ii)%species_orders%d3
      end if
      ! species list of ii atom in mlff model
      allocate(descriptor_params(ii)%species_lst(n_species),STAT=stat)
      if(stat/=0) &
        call cq_abort('acsf2b descriptor_params(ii)%species_lst: error allocating memory to n_species')
    end do
    descriptor_params%n_species = n_species

    !! Allcoate number of pairs or triplets for descriptor_params
    call allocate_acsf2b_param(descriptor_params)

    do ii=1, n_species
      ! read the two body part
      ! A-A,A-B,A-C, ...
      ! 30,30,30, ...
      read(file_id,'(a)', end=999) comment

      if (inode == ionode) then
        write(*,*) 'read_acsf2b: i, comment ', ii, comment
      end if
      ! number of pair types is equal to n_species
      n_2b_pair = descriptor_params(ii)%n_species
      allocate(dims_2b(n_2b_pair),STAT=stat)
      if(stat/=0) &
        call cq_abort('dims_2b : error allocating memory to n_2b_pair')

      ! element name
      read(file_id,'(a)') comment
      if (inode == ionode) then
        write(*,*) 'read_acsf2b: element name, comment ', comment
      end if
      ! element order
      read(file_id,'(a)') comment
      call split_string(comment, descriptor_params(ii)%species_lst, delimiter, n_species)
      if (inode == ionode) then
        write(*,'(a,a30)') 'read_acsf2b: element order, comment ', comment
        write(*,*) (descriptor_params(ii)%species_lst(jj), jj=1,n_species)
      end if

      ! Read dimensions of each pair type
      ! pair information
      !read(file_id,'(a)') comment
      ! pair dimensions
      read(file_id,*) dims_2b

      ! num_g2 = sum(nums_g2()), nums_g2_acc is index of each g2
      do jj = 1, n_2b_pair
        descriptor_params(ii)%nums_g2(jj) = dims_2b(jj)
        !! TODO: from comment to get the index of params_2b(i)
        call allocate_b2_param(descriptor_params(ii)%params_g2(jj), dims_2b(jj))
        call get_b2_param(file_id,descriptor_params(ii)%params_g2(jj), dims_2b(jj))
      end do !jj1
      if (inode == ionode) then
        write(*,*) 'read_acsf2b: after read g2 , dims_2b, n_2b_pair', dims_2b,n_2b_pair
      end if

      descriptor_params(ii)%nums_g2_acc(1) = 1
      descriptor_params(ii)%num_g2 = descriptor_params(ii)%nums_g2(1)
      if (n_2b_pair .gt. 1 ) then
        do jj = 2, n_2b_pair
          dim_2b = descriptor_params(ii)%nums_g2(jj-1)
          descriptor_params(ii)%nums_g2_acc(jj) = descriptor_params(ii)%nums_g2_acc(jj-1) + dim_2b
          dim_2b = descriptor_params(ii)%nums_g2(jj)
          descriptor_params(ii)%num_g2 = descriptor_params(ii)%num_g2 + dim_2b
        end do !jj2
      end if

      if (inode == ionode) then
        write(*,*) 'read_acsf2b: after acc_2b, n_2b_pair', descriptor_params(ii)%nums_g2_acc,n_2b_pair
      end if

      !! TODO: after read all g2,g3,g4,g5, collect all of coefficients
      num_g2 = descriptor_params(ii)%num_g2
      descriptor_params(ii)%dim_coef = num_g2
      allocate(descriptor_params(ii)%coef(descriptor_params(ii)%dim_coef),STAT=stat)
      if(stat/=0) then
        call cq_abort('acsf2b coef in descriptor_params : error allocating memory to dim_coef')
      endif
      if (inode == ionode) then
        write(*,*) 'read_acsf2b: after descriptor_params dim_coef=', descriptor_params(ii)%dim_coef
      end if

      !! collect coefficients
      shift=0
      do jj = 1, n_2b_pair
        a = shift+descriptor_params(ii)%nums_g2_acc(jj)
        b = a+descriptor_params(ii)%nums_g2(jj)-1
        if (inode == ionode) then
            write(*,*) 'read_acsf2b: in do j = 1, n_2b_pair  ', a,b,n_2b_pair
        end if
        descriptor_params(ii)%coef(a:b) = descriptor_params(ii)%params_g2(jj)%coefs &
            / descriptor_params(ii)%params_g2(jj)%scales

        do kk = 1, b-a+1
          if (inode == ionode) then
            write(*,112) descriptor_params(ii)%params_g2(jj)%eta(kk),&
                      descriptor_params(ii)%params_g2(jj)%rs(kk),&
                      descriptor_params(ii)%params_g2(jj)%scales(kk),&
                      descriptor_params(ii)%params_g2(jj)%coefs(kk)
          end if
        end do ! kk

        ! check eta type
        ! in our PCCP, we use r^2/eta^2, string as r2_div_eta2
        if (eta_type=='eta_r2') then
          if (inode == ionode) then
            write(*,*) 'This is not default eta: Do something to eta in descriptor_params', eta_type
          end if
          descriptor_params(ii)%params_g2(jj)%eta = 1.0 / sqrt(descriptor_params(ii)%params_g2(jj)%eta)
        else if (eta_type=='r2_div_eta2') then
          if (inode == ionode) then
            write(*,*) 'This is default eta: Do nothing to eta in descriptor_params',eta_type
          end if
          !descriptor_params(i)%params_g2(j)%eta = 1.0*descriptor_params(i)%params_g2(j)%eta
          !call cq_abort('Error: need to prepare this part of transfer eta')
        else if (eta_type=='r2_div_2eta2') then
          if (inode == ionode) then
            write(*,*) 'This is not default eta: Do something to eta in descriptor_params',eta_type
          end if
          descriptor_params(ii)%params_g2(jj)%eta = descriptor_params(ii)%params_g2(jj)%eta * sqrt(2.0)
          call cq_abort('Error: need to prepare this part of transfer eta')
        else
          if (inode == ionode) then
            write(*,*) 'No such eta at present',eta_type
          end if
          call cq_abort('Error: no such eta at present')
        end if
      end do !jj3
      shift=shift + descriptor_params(ii)%num_g2
      if (inode == ionode) then
        write(*,*) 'read_acsf2b: after collect coefficients '
      end if

      !! check with output
      if (inode == ionode) then
        write(*,*) 'read_acsf2b: check coef, dim_coef=',descriptor_params(ii)%dim_coef
        do jj=1, descriptor_params(ii)%dim_coef
          write(*,*) jj, descriptor_params(ii)%coef(jj)
        end do !jj4
        write(*,*) 'read_acsf2b: after check with output '
      end if

      deallocate(dims_2b)
    end do !ii1

    read(file_id,*, end=999) comment
    close(file_id)

999 if (inode == ionode) then
      print *, 'Done reading acsf2b file'
    end if
    close(file_id)

    call my_barrier
    return
112 format(4E16.8)
  end subroutine read_acsf2b
!!***

!!****f* matrix_module/allocate_split_param *
  subroutine allocate_acsf2b_param(descriptor_params)
    implicit none

    ! Passed variables
    type(acsf_param), dimension(:) :: descriptor_params

    ! Local variables
    integer :: i, n_2b_pair, dims_3b, n_spieces,stat

    n_spieces = 0
    if (descriptor_params(1)%n_species .LE. 0) then
      call cq_abort('split(1) n_species: not defined or less than 1')
    else
      n_spieces = descriptor_params(1)%n_species
    end if

    do i=1, n_spieces
      n_2b_pair = descriptor_params(i)%n_species
      !dims_3b = dims_2b * (dims_2b + 1)/2
      allocate(descriptor_params(i)%nums_g2(n_2b_pair),STAT=stat)
      if(stat/=0) &
        call cq_abort('acsf nums_g2: error allocating memory to n_2b_pair')
      allocate(descriptor_params(i)%nums_g2_acc(n_2b_pair),STAT=stat)
      if(stat/=0) &
        call cq_abort('acsf nums_g2_acc: error allocating memory to n_2b_pair')

      allocate(descriptor_params(i)%params_g2(n_2b_pair),STAT=stat)
      if(stat/=0) &
        call cq_abort('acsf params_g2: error allocating memory to n_2b_pair')
        !allocate(split(i)%nums_3b(dims_3b),STAT=stat)
        !if(stat/=0) then
        !    call cq_abort('split(1) nums_3b: error allocating memory to dims_3b')
        !endif
        !allocate(split(i)%nums_3b_acc(dims_3b),STAT=stat)
        !if(stat/=0) then
        !    call cq_abort('split(1) nums_3b_acc: error allocating memory to dims_3b')
        !endif
    end do
    return
  end subroutine allocate_acsf2b_param
!!***

!!****f* matrix_module/deallocate_split_param *
  subroutine deallocate_acsf2b_param(descriptor_params)
    implicit none

    ! Passed variables
    type(acsf_param), dimension(:), allocatable :: descriptor_params

    ! Local variables
    integer :: ii, n_2b_pair, dims_3b, n_spieces,stat

    n_spieces = 0
    if (descriptor_params(1)%n_species .LE. 0) then
      call cq_abort('split(1) n_species: not defined or less than 1')
    else
      n_spieces = descriptor_params(1)%n_species
    end if
    call start_timer(tmr_std_allocation)
    do ii=n_spieces, 1, -1
      n_2b_pair = descriptor_params(ii)%n_species
      !dims_3b = dims_2b * (dims_2b + 1)/2
      deallocate(descriptor_params(ii)%nums_g2,STAT=stat)
      if(stat/=0) &
        call cq_abort('acsf deallocate: error deallocating memory to nums_g2')
      deallocate(descriptor_params(ii)%nums_g2_acc,STAT=stat)
      if(stat/=0) &
        call cq_abort('acsf deallocate: error deallocating memory to nums_g2_acc')

      call deallocate_b2_param(descriptor_params(ii)%params_g2, n_2b_pair)

      deallocate(descriptor_params(ii)%species_lst,STAT=stat)
      if(stat/=0) &
        call cq_abort('acsf deallocate: error deallocating memory to species_lst')
    end do
    call stop_timer(tmr_std_allocation)
    return
  end subroutine deallocate_acsf2b_param
!!***


!!****f* mlff_type_module/get_feature_acsf2b *
!!
!!  NAME
!!   get_feature_acsf2b
!!  USAGE
!!
!!  PURPOSE
!!   Calculate value of descriptor for each atom
!!  INPUTS
!!   type(primary_set) :: prim
!!   type(cover_set) :: gcs
!!   type(matrix), dimension(:) :: amat
!!   real(double) :: rcut
!!   type(acsf_param), dimension(:) :: descriptor_params
!!  USES
!!
!!  AUTHOR
!!   Jianbo Lin
!!  CREATION DATE
!!   2022/07/20
!!  MODIFICATION HISTORY
!!   2022/07/27 J.Lin
!!    Modified to get pair information for machine learing
!!   2023/08/24 J.Lin
!!   Available different feature dimensions for each type of atom
!!  SOURCE
!!
  subroutine get_feature_acsf2b(prim,gcs,amat,amat_features_ML,rcut,descriptor_params)

    ! Module usage
    use datatypes
    use basic_types
    use matrix_module
    use units
    use numbers,   ONLY: pi
    use GenComms, ONLY: cq_abort, inode, ionode, my_barrier
    use global_module, ONLY: id_glob, species_glob, sf, nlpf, paof, napf, numprocs
    use group_module, ONLY: parts
    use species_module, ONLY: nsf_species, nlpf_species, npao_species, napf_species

    implicit none

    ! Passed variables
    type(primary_set) :: prim
    type(cover_set) :: gcs
    real(double) :: rcut
    type(matrix_ML), dimension (:) :: amat
    type(features_ML), dimension(:) :: amat_features_ML
    type(acsf_param), dimension(:) :: descriptor_params

    ! Local variables
    integer :: inp, cumu_ndims, neigh_spec, ip_process,ia_glob
    integer :: nn,ii,jj,kk,np,ni, ist_j, ist_k, species_orderij, i_species, j_species
    integer :: param_start, param_end, shift_dim, param_index, fp_index, gx_index

    real(double) :: rcutsq, xij, yij, zij, xik, yik, zik, rij, rik, rjk, rjk_2
    real(double) :: eta, rcut_a
    real(double) :: tmp1, tmpx, tmpy, tmpz, frc_ij, frc_ik, frc_jk
    real(double), parameter :: tol=1.0e-8_double

    ! After we have neighbor information and descriptor information
    ! Check that prim and gcs are correctly set up
    if((.NOT.ASSOCIATED(gcs%xcover)).OR. &
        (.NOT.ASSOCIATED(prim%xprim))) then
      call cq_abort('get_naba: gcs or prim without members !')
    endif
    ! rcut bohr to angstrom
    rcut_a=rcut * BohrToAng

    if (inode== ionode .and. flag_debug_mlff) then
      write(*,*) 'We are in get_naba_ML', ' id=',inode
    end if
    call my_barrier()
    ! loop over all atom pairs (atoms in primary set, max. cover set) -
    inp=1  ! Indexes primary atoms
    amat(1)%offset=0
    amat(1)%nd_offset=0
    do nn=1,prim%groups_on_node ! Partitions in primary set
      !pair check!write(*,*) 'nn loop, inode=',inode
      if(prim%nm_nodgroup(nn).gt.0) then  ! Are there atoms ?
        !write(*,*) 'if prim atoms, inode=',inode, prim%nm_nodgroup(nn) ! success until here
        !call my_barrier()
        if (inode== ionode .and. flag_debug_mlff) then
          write(*,*) 'We are in get_naba_ML prim%nm_nodgroup(nn).gt.0'
        end if
        if (amat(nn)%n_atoms .ne. prim%nm_nodgroup(nn)) then ! Redundant, but useful
          call cq_abort('amat(nn)%n_atoms not equal prim%nm_nodgroup(nn)')
        end if
        !write(io_lun,*) 'Starting group with atoms: ',nn,prim%nm_nodgroup(nn)
        do ii=1,prim%nm_nodgroup(nn)  ! Loop over atoms in partition
          !write(*,*) 'j loop, inode=',inode
          !call my_barrier()
          !Species of i is  amat(nn)%i_species(i)
          i_species = amat(nn)%i_species(ii)
          ia_glob=prim%ig_prim(prim%nm_nodbeg(nn)+ii-1)
          do jj=1,  amat(nn)%n_nab(ii)
            ist_j = amat(nn)%i_acc(ii)+jj-1
            j_species = amat(nn)%j_species(ist_j)
            species_orderij = descriptor_params(i_species)%species_orders%d2(j_species)

            rij=amat(nn)%radius(ist_j) * BohrToAng
            xij=amat(nn)%dx(ist_j) * BohrToAng
            yij=amat(nn)%dy(ist_j) * BohrToAng
            zij=amat(nn)%dz(ist_j) * BohrToAng

            !! function of cut off
            frc_ij = 0.5 * (cos(pi * rij / rcut_a) + 1)

            !! Start calculate the two body terms
            shift_dim = 0
            param_start = descriptor_params(i_species)%nums_g2_acc(species_orderij)
            param_end = param_start + descriptor_params(i_species)%nums_g2(species_orderij) - 1
            !!check params
            if (inode== ionode .and. ii==1 .and. flag_debug_mlff) then
              write(*,102) param_start, param_end, i_species,j_species,amat(nn)%n_nab(ii)
            end if
            do param_index = param_start, param_end
              fp_index = param_index + shift_dim
              ! Todo : check if species_orderij is correct
              gx_index = param_index - param_start + 1
              eta = descriptor_params(i_species)%params_g2(species_orderij)%eta(gx_index)

              tmp1 = frc_ij * exp(- (rij/eta) ** 2) / rij
              amat_features_ML(nn)%id_atom(ii)%fpx(fp_index) = &
                  amat_features_ML(nn)%id_atom(ii)%fpx(fp_index) + xij*tmp1
              amat_features_ML(nn)%id_atom(ii)%fpy(fp_index) = &
                  amat_features_ML(nn)%id_atom(ii)%fpy(fp_index) + yij*tmp1
              amat_features_ML(nn)%id_atom(ii)%fpz(fp_index) = &
                  amat_features_ML(nn)%id_atom(ii)%fpz(fp_index) + zij*tmp1
              !if (inode== ionode .and. ia_glob==1) then
              !    write(1986,101)  eta, fp_index,  tmp1, xij, xij*tmp1, amat_features_ML(nn)%fpx(i, fp_index)
              !end if
            end do ! two-body terms
            !pair check!write(*,*) 'after part_nd_nabs in get_feature_acsf2b', inode
          end do ! j, two-body

          inp=inp+1  ! Indexes primary-set atoms
        enddo ! End prim%nm_nodgroup
        !pair check!write(*,*) 'after prim%nm_nodgroup', inode

      else
        if(flag_debug_mlff) &
            write(*, *) 'Warning: No atoms in this partition get_feature_acsf', inode, nn
      endif ! End if(prim%nm_nodgroup>0)
    enddo ! End part_on_node
    100 format(6i8,4e25.16)
    101 format(e25.16, i8,4e25.16)
    102 format('check fp index,start,end',2i8,'species=',3i8)
  end subroutine get_feature_acsf2b
!!***


!!****f* matrix_module/allocate_split_param *
  subroutine allocate_split_param(descriptor_params)
    implicit none

    ! Passed variables
    type(split_param), dimension(:) :: descriptor_params

    ! Local variables
    integer :: i, n_2b, n_3b, n_spieces,stat

    n_spieces = 0
    if (descriptor_params(1)%n_species .LE. 0) then
      call cq_abort('split(1) n_species: not defined or less than 1')
    else
      n_spieces = descriptor_params(1)%n_species
    end if

    do i=1, n_spieces
      ! Two body part
      ! number of pair types
      n_2b = descriptor_params(i)%n_species
      allocate(descriptor_params(i)%nums_2b(n_2b),STAT=stat)
      if(stat/=0) &
        call cq_abort('split nums_2b: error allocating memory to n_2b')
      allocate(descriptor_params(i)%nums_2b_acc(n_2b),STAT=stat)
      if(stat/=0) &
        call cq_abort('split nums_2b_acc: error allocating memory to n_2b')
      allocate(descriptor_params(i)%params_2b(n_2b),STAT=stat)
      if(stat/=0) &
        call cq_abort('split params_2b: error allocating memory to n_2b')

      ! Three body part
      ! number of triplet types
      n_3b = n_2b * (n_2b + 1)/2
      allocate(descriptor_params(i)%nums_3b(n_3b),STAT=stat)
      if(stat/=0) &
        call cq_abort('split nums_3b: error allocating memory to n_3b')
      allocate(descriptor_params(i)%nums_3b_acc(n_3b),STAT=stat)
      if(stat/=0) &
        call cq_abort('split nums_3b_acc: error allocating memory to n_3b')
      allocate(descriptor_params(i)%params_3b(n_3b),STAT=stat)
      if(stat/=0) &
          call cq_abort('split params_3b: error allocating memory to n_3b')
    end do
    return
  end subroutine allocate_split_param
!!***

!!****f* matrix_module/allocate_split_param *
  subroutine deallocate_split_param(descriptor_params)
    implicit none

    ! Passed variables
    type(split_param), dimension(:), allocatable :: descriptor_params

    ! Local variables
    integer :: i, n_2b, n_3b, n_spieces,stat

    n_spieces = 0
    if (descriptor_params(1)%n_species .LE. 0) then
        call cq_abort('split(1) n_species: not defined or less than 1')
    else
        n_spieces = descriptor_params(1)%n_species
    end if
    call start_timer(tmr_std_allocation)
    do i=n_spieces, 1,-1
      ! Two body part
      ! number of pair types
      n_2b = descriptor_params(i)%n_species
      deallocate(descriptor_params(i)%nums_2b,STAT=stat)
      if(stat/=0) &
        call cq_abort('split deallocate: error deallocating memory to nums_2b')
      deallocate(descriptor_params(i)%nums_2b_acc,STAT=stat)
      if(stat/=0) &
        call cq_abort('split deallocate: error deallocating memory to nums_2b_acc')

      call deallocate_b2_param(descriptor_params(i)%params_2b,n_2b)

      ! Three body part
      ! number of triplet types
      n_3b = n_2b * (n_2b + 1)/2
      deallocate(descriptor_params(i)%nums_3b,STAT=stat)
      if(stat/=0) &
        call cq_abort('split deallocate: error deallocating memory to nums_3b')
      deallocate(descriptor_params(i)%nums_3b_acc,STAT=stat)
      if(stat/=0) &
        call cq_abort('split deallocate: error deallocating memory to nums_3b_acc')

      call deallocate_b3_param(descriptor_params(i)%params_3b,n_3b)
    end do
    call stop_timer(tmr_std_allocation)
    return
  end subroutine deallocate_split_param
!!***

!!****f*  ML_type/read_split *
!!
!!  NAME
!!   read_split
!!  PURPOSE
!!   read model information of machine learning for CONQUEST
!!  AUTHOR
!!   Jianbo Lin
!!  CREATION DATE
!!   2022/03/14
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine read_split(filename, descriptor_params)

    use GenComms, ONLY: cq_abort, inode, ionode,my_barrier
    use global_module, ONLY: area_integn
    use memory_module, ONLY: reg_alloc_mem, type_int

    implicit none

    ! Passed variables
    type(split_param),allocatable  :: descriptor_params(:)
    character(len=*)           :: filename

    ! Local variables
    integer                    :: ii, jj, kk, n_species, n_2b, n_3b, dim_2b, dim_3b
    integer                    :: a,b,shift, num_2b, num_3b
    integer                    :: stat
    integer                    :: file_id=2022
    integer, allocatable       :: dims_2b(:), dims_3b(:)
    character(len=180)         :: comment
    character(len=30)          :: eta_type, delimiter

    delimiter=' '
    open(unit=file_id,file=filename,iostat=stat,status='old')
    read(file_id,*, end=999) comment
    read(file_id,*, end=999) comment
    read(file_id,*, end=999) comment
    read(file_id,*) n_species
    !! TODO: remove in the future about eta type
    read(file_id,*, end=999) comment
    read(file_id,*) eta_type
    eta_type=TRIM(eta_type)

    !! Allocate n_spieces of descriptor_params
    allocate(descriptor_params(n_species),STAT=stat)
    if(stat/=0) then
      call cq_abort('split2b3b descriptor_params: error allocating memory to n_species')
    endif
    !! TODO: consider a well-designed species order
    do ii=1, n_species
      ! species_orders only depending on input of species order
      !call ini_species_order(descriptor_params(ii)%species_orders, n_species)
      ! species_orders depending on center atom according to the input of species order
      call ini_species_order_center(descriptor_params(ii)%species_orders, n_species, ii)
      if (inode == ionode) then
        write(*,*) 'read_split2b3b: species_order d2 ', descriptor_params(ii)%species_orders%d2
        write(*,*) 'read_split2b3b: species_order d3', descriptor_params(ii)%species_orders%d3
      end if
      ! species list of ii atom in mlff model
      allocate(descriptor_params(ii)%species_lst(n_species),STAT=stat)
    end do ! ii1
    descriptor_params%n_species = n_species

    !! Allcoate number of pairs or triplets for descriptor_params
    call allocate_split_param(descriptor_params)

    do ii=1, n_species
      ! read the two body part
      ! A-A,A-B,A-C, ...
      ! 30,30,30, ...
      ! read comment lines
      read(file_id,'(a)', end=999) comment !To seperate from other parameters and show detail of descriptor
      read(file_id,'(a)', end=999) comment !To seperate from other parameters and show detail of descriptor

      if (inode == ionode) then
          write(*,*) 'read_split2b3b: ii, comment ', ii, comment
      end if

      ! read center atom, and order of species
      ! element name
      read(file_id,'(a)') comment
      if (inode == ionode) then
        write(*,*) 'read_split2b3b: element name, comment ', comment
      end if
      ! element order
      read(file_id,'(a)') comment
      call split_string(comment, descriptor_params(ii)%species_lst, delimiter, n_species)
      if (inode == ionode) then
        write(*,'(a,a30)') 'read_split2b3b: element order, comment ', comment
        write(*,*) (descriptor_params(ii)%species_lst(jj), jj=1,n_species)
      end if

      ! Read dimensions of each pair type
      ! number of pair types is equal to n_species
      n_2b = descriptor_params(ii)%n_species
      allocate(dims_2b(n_2b),STAT=stat)
      if(stat/=0) &
        call cq_abort('dims_2b : error allocating memory to n_2b')

      n_3b = n_2b * (n_2b+1)/2
      allocate(dims_3b(n_3b),STAT=stat)
      if(stat/=0) &
        call cq_abort('dims_3b : error allocating memory to n_3b')

      !! Read dimensions of two body terms
      ! pair information
      !read(file_id,'(a)') comment
      read(file_id,*) dims_2b
      !! Read dimensions of three body terms
      ! triplet information
      !read(file_id,'(a)') comment
      read(file_id,*) dims_3b

      !! Read two body terms
      ! num_g2 = sum(nums_g2()), nums_g2_acc is index of each g2
      do jj = 1, n_2b
        descriptor_params(ii)%nums_2b(jj) = dims_2b(jj)
        !! TODO: from comment to get the index of params_2b(i)
        call allocate_b2_param(descriptor_params(ii)%params_2b(jj), dims_2b(jj))
        call get_b2_param(file_id,descriptor_params(ii)%params_2b(jj), dims_2b(jj))
      end do !jj1
      if (inode == ionode) then
        write(*,*) 'read_split2b3b: after read 2b , dims_2b, n_2b ', dims_2b,n_2b
      end if

      !! update index of two body terms
      descriptor_params(ii)%nums_2b_acc(1) = 1
      descriptor_params(ii)%num_2b = descriptor_params(ii)%nums_2b(1)
      if (n_2b .gt. 1 ) then
        do jj = 2, n_2b
          dim_2b = descriptor_params(ii)%nums_2b(jj-1)
          descriptor_params(ii)%nums_2b_acc(jj) = descriptor_params(ii)%nums_2b_acc(jj-1) + dim_2b
          dim_2b = descriptor_params(ii)%nums_2b(jj)
          descriptor_params(ii)%num_2b = descriptor_params(ii)%num_2b + dim_2b
        end do !jj2
      end if

      if (inode == ionode) then
        write(*,*) 'read_split2b3b: after acc_2b, n2b ', descriptor_params(ii)%nums_2b_acc,n_2b
      end if

      !! Read three body terms
      do jj = 1, n_3b
        descriptor_params(ii)%nums_3b(jj) = dims_3b(jj)
        !! TODO: from comment to get the index of params_2b(ii)
        call allocate_b3_param(descriptor_params(ii)%params_3b(jj), dims_3b(jj))
        call get_b3_param(file_id,descriptor_params(ii)%params_3b(jj), dims_3b(jj))
      end do !jj3,n_3b
      if (inode == ionode) then
        write(*,*) 'read_split2b3b: after read 3b , dims_3b, n_3b ', dims_3b,n_3b
      end if

      !! update index of three body terms
      descriptor_params(ii)%nums_3b_acc(1) = 1
      descriptor_params(ii)%num_3b = descriptor_params(ii)%nums_3b(1)
      if (n_3b .gt. 1 ) then
        do jj = 2, n_3b
          dim_3b = descriptor_params(ii)%nums_3b(jj-1)
          descriptor_params(ii)%nums_3b_acc(jj) = descriptor_params(ii)%nums_3b_acc(jj-1) + dim_3b
          dim_3b = descriptor_params(ii)%nums_3b(jj)
          descriptor_params(ii)%num_3b = descriptor_params(ii)%num_3b + dim_3b
        end do !jj4,n_3b
      end if

      if (inode == ionode) then
        write(*,*) 'read_split2b3b: after acc, n_3b ', descriptor_params(ii)%nums_3b_acc,n_3b
      end if
      !###############

      !! after read all two body terms and three body terms, collect all of coefficients
      num_2b = descriptor_params(ii)%num_2b
      num_3b = descriptor_params(ii)%num_3b
      descriptor_params(ii)%dim_coef = num_2b + num_3b
      allocate(descriptor_params(ii)%coef(descriptor_params(ii)%dim_coef),STAT=stat)
      if(stat/=0) &
        call cq_abort('split2b3b coef in descriptor_params : error allocating memory to dim_coef')
      if (inode == ionode) then
        write(*,*) 'read_split2b3b: after descriptor_params dim_coef= ', descriptor_params(ii)%dim_coef
      end if

      !! Collect coefficients
      !! collect coefficients of two body terms
      shift=0
      do jj = 1, n_2b
        a = shift+descriptor_params(ii)%nums_2b_acc(jj)
        b = a+descriptor_params(ii)%nums_2b(jj)-1
        if (inode == ionode) then
          write(*,*) 'read_split2b3b: in do j = 1, n_2b  ', a,b,n_2b
        end if
        descriptor_params(ii)%coef(a:b) = descriptor_params(ii)%params_2b(jj)%coefs &
            / descriptor_params(ii)%params_2b(jj)%scales

        !check parameters of three body terms
        do kk = 1, b-a+1
          if (inode == ionode) then
            write(*,114) descriptor_params(ii)%params_2b(jj)%eta(kk),&
                      descriptor_params(ii)%params_2b(jj)%rs(kk),&
                      descriptor_params(ii)%params_2b(jj)%scales(kk),&
                      descriptor_params(ii)%params_2b(jj)%coefs(kk)
          end if
        end do !kk1

        ! check eta type, and transfter to standard form
        ! in our PCCP, we use r^2/eta^2, string as r2_div_eta2
        if (eta_type=='eta_r2') then
          if (inode == ionode) &
            write(*,*) 'This is not default eta: Do something to eta in descriptor_params', eta_type
          descriptor_params(ii)%params_2b(jj)%eta = 1/ sqrt(descriptor_params(ii)%params_2b(jj)%eta)
        else if (eta_type=='r2_div_eta2') then
          if (inode == ionode) &
            write(*,*) 'This is default eta: Do nothing to eta in descriptor_params',eta_type
        else if (eta_type=='r2_div_2eta2') then
          if (inode == ionode) then
            write(*,*) 'This is not default eta: Do something to eta in descriptor_params',eta_type
          end if
          descriptor_params(ii)%params_2b(jj)%eta = descriptor_params(ii)%params_2b(jj)%eta * sqrt(2.0)
          call cq_abort('Error: need to prepare this part of transfer eta')
        else
          if (inode == ionode) then
            write(*,*) 'No such eta at present',eta_type
          end if
          call cq_abort('Error: no such eta at present')
        end if
      end do !jj5 n_2b
      !! collect coefficients of three body terms
      shift=shift + descriptor_params(ii)%num_2b
      do jj = 1, n_3b
        a = shift+descriptor_params(ii)%nums_3b_acc(jj)
        b = a+descriptor_params(ii)%nums_3b(jj)-1
        descriptor_params(ii)%coef(a:b) = descriptor_params(ii)%params_3b(jj)%coefs &
            / descriptor_params(ii)%params_3b(jj)%scales

        !check parameters of three body terms
        if (inode == ionode) then
          write(*,*) 'read_split2b3b: in do j = 1, n_3b  ', a,b,n_3b
          do kk = 1, b-a+1
            write(*,115) descriptor_params(ii)%params_3b(jj)%eta1(kk),&
                      descriptor_params(ii)%params_3b(jj)%eta2(kk),&
                      descriptor_params(ii)%params_3b(jj)%eta3(kk),&
                      descriptor_params(ii)%params_3b(jj)%scales(kk),&
                      descriptor_params(ii)%params_3b(jj)%coefs(kk)
          end do ! kk
        end if

        ! check eta type, and transfter to standard form
        ! in our PCCP, we use r^2/eta^2, string as r2_div_eta2
        if (eta_type=='eta_r2') then
          if (inode == ionode) then
            write(*,*) 'This is not default eta: Do something to eta in descriptor_params', eta_type
          end if
          descriptor_params(ii)%params_3b(jj)%eta1 = 1/ sqrt(descriptor_params(ii)%params_3b(jj)%eta1)
          descriptor_params(ii)%params_3b(jj)%eta2 = 1/ sqrt(descriptor_params(ii)%params_3b(jj)%eta2)
          descriptor_params(ii)%params_3b(jj)%eta3 = 1/ sqrt(descriptor_params(ii)%params_3b(jj)%eta3)
        else if (eta_type=='r2_div_eta2') then
          if (inode == ionode) &
            write(*,*) 'This is default eta: Do nothing to eta in descriptor_params',eta_type
        else if (eta_type=='r2_div_2eta2') then
          if (inode == ionode) then
            write(*,*) 'This is not default eta: Do something to eta in descriptor_params',eta_type
          end if
          descriptor_params(ii)%params_3b(jj)%eta1 = sqrt(2.0) * descriptor_params(ii)%params_3b(jj)%eta1
          descriptor_params(ii)%params_3b(jj)%eta2 = sqrt(2.0) * descriptor_params(ii)%params_3b(jj)%eta2
          descriptor_params(ii)%params_3b(jj)%eta3 = sqrt(2.0) * descriptor_params(ii)%params_3b(jj)%eta3
        else
          if (inode == ionode) then
            write(*,*) 'No such eta at present',eta_type
          end if
          call cq_abort('Error: no such eta type at present')
        end if
      end do !jj6 n_3b

      !! check with output
      if (inode == ionode) then
        write(*,*) 'read_split2b3b: after collect coefficients, dim_coef=',descriptor_params(ii)%dim_coef
        write(*,*) 'read_split2b3b: start check with output '
        do jj=1, descriptor_params(ii)%dim_coef
          write(*,*) jj, descriptor_params(ii)%coef(jj)
        end do !jj7
        write(*,*) 'read_split2b3b: after check with output '
      end if

      deallocate(dims_2b)
      deallocate(dims_3b)
    end do !ii n_species

    read(file_id,*, end=999) comment
    close(file_id)

999 if (inode == ionode) then
      print *, 'Done reading split2b3b file'
    end if
    close(file_id)

    call my_barrier
    return
114 format(4E16.8)
115 format(5E16.8)
  end subroutine read_split
!!***



!!****f* mlff_type_module/get_feature_ML *
!!
!!  NAME
!!   get_feature_ML
!!  USAGE
!!
!!  PURPOSE
!!   Creates neighbour lists for a matrix range rcut
!!  INPUTS
!!   type(primary_set) :: prim
!!   type(cover_set) :: gcs
!!   type(matrix), dimension(:) :: amat
!!   real(double) :: rcut
!!  USES
!!
!!  AUTHOR
!!   Jianbo Lin
!!  CREATION DATE
!!   2022/07/20
!!  MODIFICATION HISTORY
!!   2022/07/27 J.Lin
!!    Modified to get pair information for machine learing
!!   2023/08/24 J.Lin
!!   Available different feature dimensions for each type of atom
!!   2023/08/27 J.Lin
!!   Correction for multi-element bond order
!!  SOURCE
!!
  subroutine get_feature_split(prim,gcs,amat,amat_features_ML,rcut,descriptor_params)

    ! Module usage
    use datatypes

    use basic_types
    use matrix_module
    use units
    use numbers,   ONLY: pi,very_small
    use GenComms, ONLY: cq_abort, inode, ionode, my_barrier
    use global_module, ONLY: id_glob, species_glob, sf, nlpf, paof, napf, numprocs
    use group_module, ONLY: parts
    use species_module, ONLY: nsf_species, nlpf_species, npao_species, napf_species

    implicit none

    ! Passed variables
    type(primary_set) :: prim
    type(cover_set) :: gcs
    real(double) :: rcut
    type(matrix_ML), dimension (:) :: amat
    type(features_ML), dimension(:) :: amat_features_ML
    type(split_param), dimension(:) :: descriptor_params

    ! Local variables
    integer :: inp, cumu_ndims, neigh_spec, ip_process
    integer :: nn,ii,jj,kk, ist_j, ist_k, ist_ijk, i_3b, ia_glob
    integer :: species_order2b, species_order3b, i_species, j_species, k_species, tmp_species
    integer :: param_start, param_end, shift_dim, param_index, fp_index, gx_index

    real(double) :: xij0, yij0, zij0, xik0, yik0, zik0, rij0, rik0, rjk0 ! from ML matrix
    !real(double) :: xij, yij, zij, xik, yik, zik, rij, rik, rjk          ! special bond order
    real(double) :: eta, rs, eta1, eta2, eta3, rcut_2, rcut_a, rcutsq
    real(double) :: tmp1, tmpx, tmpy, tmpz, tmpr, frc_ij, frc_ik, frc_jk, frc_3b
    real(double) :: tmp123, tmp132, tmp213, tmp231, tmp312, tmp321
    real(double) :: proj_x,proj_y,proj_z
    real(double), parameter :: tol=1.0e-8_double
    real(double) :: feature_1 ! check feature

    ! After we have neighbor information and descriptor information
    ! Check that prim and gcs are correctly set up
    if((.NOT.ASSOCIATED(gcs%xcover)).OR. &
        (.NOT.ASSOCIATED(prim%xprim))) then
      call cq_abort('get_naba: gcs or prim without members !')
    endif
    ! In this subrputine, Units: Length=Angstrom, Energy=eV, Force=eV/Anstrom
    ! rcut bohr to angstrom
    rcut_a=rcut * BohrToAng
    rcut_2=rcut_a * rcut_a

    call my_barrier()
    ! loop over all atom pairs (atoms in primary set, max. cover set) -
    inp=1  ! Indexes primary atoms
    do nn=1,prim%groups_on_node ! Partitions in primary set
      !pair check!write(*,*) 'nn loop, inode=',inode
      if(prim%nm_nodgroup(nn).gt.0) then  ! Are there atoms ?
        !write(*,*) 'if prim atoms, inode=',inode, prim%nm_nodgroup(nn) ! success until here
        !call my_barrier()
        if (amat(nn)%n_atoms .ne. prim%nm_nodgroup(nn)) then ! Redundant, but useful
          call cq_abort('amat(nn)%n_atoms not equal prim%nm_nodgroup(nn)')
        end if
        !write(io_lun,*) 'Starting group with atoms: ',nn,prim%nm_nodgroup(nn)
        do ii=1,prim%nm_nodgroup(nn)  ! Loop over atoms in partition
          ! From nn and ii, glob id and species of i-th atom are shown as following
          !Species of i is  amat(nn)%i_species(i)
          i_species = amat(nn)%i_species(ii)
          ia_glob = prim%ig_prim(prim%nm_nodbeg(nn)+ii-1)
          do jj=1,  amat(nn)%n_nab(ii)
            ist_j = amat(nn)%i_acc(ii)+jj-1
            j_species = amat(nn)%j_species(ist_j)
            species_order2b = descriptor_params(i_species)%species_orders%d2(j_species)

            rij0=amat(nn)%radius(ist_j) * BohrToAng
            xij0=amat(nn)%dx(ist_j) * BohrToAng
            yij0=amat(nn)%dy(ist_j) * BohrToAng
            zij0=amat(nn)%dz(ist_j) * BohrToAng

            !! function of cut off
            frc_ij = 0.5 * (cos(pi * rij0 / rcut_a) + 1)

            !! Start calculate the two body terms
            shift_dim = 0
            param_start = descriptor_params(i_species)%nums_2b_acc(species_order2b)
            param_end = param_start + descriptor_params(i_species)%nums_2b(species_order2b) - 1
            !!check params
            if (ia_glob==1 .and. flag_debug_mlff) then
              write(*,102) param_start, param_end, i_species,j_species, species_order2b, amat(nn)%n_nab(ii)
            end if

            do param_index = param_start, param_end
              fp_index = param_index + shift_dim
              ! Todo: check dimension of params_2b(j)
              gx_index = param_index - param_start + 1
              eta = descriptor_params(i_species)%params_2b(species_order2b)%eta(gx_index)
              eta = 1.0/eta**2
              rs = rij0 - descriptor_params(i_species)%params_2b(species_order2b)%rs(gx_index)

              tmp1 = frc_ij * exp(- eta * rs ** 2 )  !* eta * rs / rij
              amat_features_ML(nn)%id_atom(ii)%fpx(fp_index) = &
                  amat_features_ML(nn)%id_atom(ii)%fpx(fp_index) + xij0 * tmp1
              amat_features_ML(nn)%id_atom(ii)%fpy(fp_index) = &
                  amat_features_ML(nn)%id_atom(ii)%fpy(fp_index) + yij0 * tmp1
              amat_features_ML(nn)%id_atom(ii)%fpz(fp_index) = &
                  amat_features_ML(nn)%id_atom(ii)%fpz(fp_index) + zij0 * tmp1
            end do ! two-body terms

            do kk=1,  amat(nn)%n_triplet(ist_j)
              ist_ijk=amat(nn)%triplet_acc(ist_j)+kk-1
              ist_k=amat(nn)%ist_triplet(ist_ijk)
              k_species = amat(nn)%j_species(ist_k)

              rik0=amat(nn)%radius(ist_k) * BohrToAng
              xik0=amat(nn)%dx(ist_k) * BohrToAng
              yik0=amat(nn)%dy(ist_k) * BohrToAng
              zik0=amat(nn)%dz(ist_k) * BohrToAng

              !! Todo: check order for ijk
              species_order3b = descriptor_params(i_species)%species_orders%d3(j_species, k_species)
              rjk0 = (xij0-xik0)**2+(yij0-yik0)**2+(zij0-zik0)**2
              rjk0 = sqrt(rjk0)

              !! function of cut off
              frc_ik = 0.5 * (cos(pi * rik0 / rcut_a) + 1)
              frc_jk = 0.5 * (cos(pi * rjk0 / rcut_a) + 1)
              frc_3b = frc_ij * frc_ik * frc_jk

              !! Start calculate the three body terms
              shift_dim = descriptor_params(i_species)%num_2b
              param_start = descriptor_params(i_species)%nums_3b_acc(species_order3b)
              param_end = param_start + descriptor_params(i_species)%nums_3b(species_order3b) - 1

              if (ia_glob==1 .and. flag_debug_mlff) then
                write(*,103) jj, kk, shift_dim,param_start, param_end, &
                    i_species,j_species,k_species, species_order3b, amat(nn)%n_nab(ii)
                feature_1 = amat_features_ML(nn)%id_atom(ii)%fpx(121)
              end if

              ! Checking X-XX, X-AA, X-XA, X-AB
              if (j_species == k_species) then
                ! X-XX or X-AA
                if (i_species == j_species) then
                  ! X-XX
                  if (ia_glob==1 .and. flag_debug_mlff) then
                    write(*,104) param_start, param_end, &
                      i_species,j_species,k_species,amat(nn)%n_nab(ii), 'X-XX'
                  end if
                  ! debug
                  !cycle
                  do param_index= param_start, param_end
                    fp_index = param_index + shift_dim

                    ! Todo: check dimension of params_2b(j)
                    gx_index = param_index - param_start + 1
                    eta1 = descriptor_params(i_species)%params_3b(species_order3b)%eta1(gx_index)
                    eta2 = descriptor_params(i_species)%params_3b(species_order3b)%eta2(gx_index)
                    eta3 = descriptor_params(i_species)%params_3b(species_order3b)%eta3(gx_index)

                    ! transfter for efficiency
                    eta1 = 1.0/eta1 ** 2
                    eta2 = 1.0/eta2 ** 2
                    eta3 = 1.0/eta3 ** 2

                    ! Terms for A-A-A
                    tmp123 = exp(-eta1 * rij0 ** 2 - eta2 * rik0 ** 2 - eta3 * rjk0 ** 2)
                    tmp132 = exp(-eta1 * rij0 ** 2 - eta3 * rik0 ** 2 - eta2 * rjk0 ** 2)
                    tmp213 = exp(-eta2 * rij0 ** 2 - eta1 * rik0 ** 2 - eta3 * rjk0 ** 2)
                    tmp231 = exp(-eta2 * rij0 ** 2 - eta3 * rik0 ** 2 - eta1 * rjk0 ** 2)
                    tmp312 = exp(-eta3 * rij0 ** 2 - eta1 * rik0 ** 2 - eta2 * rjk0 ** 2)
                    tmp321 = exp(-eta3 * rij0 ** 2 - eta2 * rik0 ** 2 - eta1 * rjk0 ** 2)

                    proj_x = 0.0
                    proj_y = 0.0
                    proj_z = 0.0
                    proj_x = proj_x + tmp123 * (xij0 * eta1 + xik0 * eta2)
                    proj_x = proj_x + tmp132 * (xij0 * eta1 + xik0 * eta3)
                    proj_x = proj_x + tmp213 * (xij0 * eta2 + xik0 * eta1)
                    proj_x = proj_x + tmp231 * (xij0 * eta2 + xik0 * eta3)
                    proj_x = proj_x + tmp312 * (xij0 * eta3 + xik0 * eta1)
                    proj_x = proj_x + tmp321 * (xij0 * eta3 + xik0 * eta2)

                    proj_y = proj_y + tmp123 * (yij0 * eta1 + yik0 * eta2)
                    proj_y = proj_y + tmp132 * (yij0 * eta1 + yik0 * eta3)
                    proj_y = proj_y + tmp213 * (yij0 * eta2 + yik0 * eta1)
                    proj_y = proj_y + tmp231 * (yij0 * eta2 + yik0 * eta3)
                    proj_y = proj_y + tmp312 * (yij0 * eta3 + yik0 * eta1)
                    proj_y = proj_y + tmp321 * (yij0 * eta3 + yik0 * eta2)

                    proj_z = proj_z + tmp123 * (zij0 * eta1 + zik0 * eta2)
                    proj_z = proj_z + tmp132 * (zij0 * eta1 + zik0 * eta3)
                    proj_z = proj_z + tmp213 * (zij0 * eta2 + zik0 * eta1)
                    proj_z = proj_z + tmp231 * (zij0 * eta2 + zik0 * eta3)
                    proj_z = proj_z + tmp312 * (zij0 * eta3 + zik0 * eta1)
                    proj_z = proj_z + tmp321 * (zij0 * eta3 + zik0 * eta2)

                    tmpx = frc_3b * proj_x
                    tmpy = frc_3b * proj_y
                    tmpz = frc_3b * proj_z

                    amat_features_ML(nn)%id_atom(ii)%fpx(fp_index) = &
                        amat_features_ML(nn)%id_atom(ii)%fpx(fp_index) + tmpx
                    amat_features_ML(nn)%id_atom(ii)%fpy(fp_index) = &
                        amat_features_ML(nn)%id_atom(ii)%fpy(fp_index) + tmpy
                    amat_features_ML(nn)%id_atom(ii)%fpz(fp_index) = &
                        amat_features_ML(nn)%id_atom(ii)%fpz(fp_index) + tmpz
                  end do ! three-body terms
                else
                  !! X-AA
                  if (ia_glob<1 .and. flag_debug_mlff) then
                  write(*,104) param_start, param_end, &
                      i_species,j_species,k_species,amat(nn)%n_nab(ii), 'X-AA'
                  end if
                  ! debug
                  !cycle
                  do param_index= param_start, param_end
                    fp_index = param_index + shift_dim

                    ! Todo: check dimension of params_2b(j)
                    gx_index = param_index - param_start + 1
                    eta1 = descriptor_params(i_species)%params_3b(species_order3b)%eta1(gx_index)
                    eta2 = descriptor_params(i_species)%params_3b(species_order3b)%eta2(gx_index)
                    eta3 = descriptor_params(i_species)%params_3b(species_order3b)%eta3(gx_index)

                    ! transfter for efficiency
                    eta1 = 1.0/eta1 ** 2
                    eta2 = 1.0/eta2 ** 2
                    eta3 = 1.0/eta3 ** 2

                    ! Terms for X-AA
                    tmp123 = exp(-eta1 * rij0 ** 2 - eta2 * rik0 ** 2 - eta3 * rjk0 ** 2)
                    tmp213 = exp(-eta2 * rij0 ** 2 - eta1 * rik0 ** 2 - eta3 * rjk0 ** 2)

                    proj_x = 0.0
                    proj_y = 0.0
                    proj_z = 0.0
                    proj_x = proj_x + tmp123 * (xij0 * eta1 + xik0 * eta2)
                    proj_x = proj_x + tmp213 * (xij0 * eta2 + xik0 * eta1)

                    proj_y = proj_y + tmp123 * (yij0 * eta1 + yik0 * eta2)
                    proj_y = proj_y + tmp213 * (yij0 * eta2 + yik0 * eta1)

                    proj_z = proj_z + tmp123 * (zij0 * eta1 + zik0 * eta2)
                    proj_z = proj_z + tmp213 * (zij0 * eta2 + zik0 * eta1)

                    tmpx = frc_3b * proj_x
                    tmpy = frc_3b * proj_y
                    tmpz = frc_3b * proj_z

                    amat_features_ML(nn)%id_atom(ii)%fpx(fp_index) = &
                        amat_features_ML(nn)%id_atom(ii)%fpx(fp_index) + tmpx
                    amat_features_ML(nn)%id_atom(ii)%fpy(fp_index) = &
                        amat_features_ML(nn)%id_atom(ii)%fpy(fp_index) + tmpy
                    amat_features_ML(nn)%id_atom(ii)%fpz(fp_index) = &
                        amat_features_ML(nn)%id_atom(ii)%fpz(fp_index) + tmpz
                  end do ! three-body terms
                end if ! X-XX or X-AA
              else
                ! X-XA, X-AX, X-AB, X-BA
                if (i_species == j_species) then
                  ! X-XA
                  if (ia_glob==1 .and. flag_debug_mlff) then
                    write(*,104) param_start, param_end, &
                      i_species,j_species,k_species,amat(nn)%n_nab(ii), 'X-XA'
                  end if

                  do param_index= param_start, param_end
                    fp_index = param_index + shift_dim

                    ! Todo: check dimension of params_2b(j)
                    gx_index = param_index - param_start + 1
                    eta1 = descriptor_params(i_species)%params_3b(species_order3b)%eta1(gx_index)
                    eta2 = descriptor_params(i_species)%params_3b(species_order3b)%eta2(gx_index)
                    eta3 = descriptor_params(i_species)%params_3b(species_order3b)%eta3(gx_index)

                    ! transfter for efficiency
                    eta1 = 1.0/eta1 ** 2
                    eta2 = 1.0/eta2 ** 2
                    eta3 = 1.0/eta3 ** 2

                    ! Terms for A-A-A
                    tmp123 = exp(-eta1 * rij0 ** 2 - eta2 * rik0 ** 2 - eta3 * rjk0 ** 2)
                    tmp132 = exp(-eta1 * rij0 ** 2 - eta3 * rik0 ** 2 - eta2 * rjk0 ** 2)

                    proj_x = 0.0
                    proj_y = 0.0
                    proj_z = 0.0
                    proj_x = proj_x + tmp123 * (xij0 * eta1 + xik0 * eta2)
                    proj_x = proj_x + tmp132 * (xij0 * eta1 + xik0 * eta3)

                    proj_y = proj_y + tmp123 * (yij0 * eta1 + yik0 * eta2)
                    proj_y = proj_y + tmp132 * (yij0 * eta1 + yik0 * eta3)

                    proj_z = proj_z + tmp123 * (zij0 * eta1 + zik0 * eta2)
                    proj_z = proj_z + tmp132 * (zij0 * eta1 + zik0 * eta3)

                    tmpx = frc_3b * proj_x
                    tmpy = frc_3b * proj_y
                    tmpz = frc_3b * proj_z

                    amat_features_ML(nn)%id_atom(ii)%fpx(fp_index) = &
                        amat_features_ML(nn)%id_atom(ii)%fpx(fp_index) + tmpx
                    amat_features_ML(nn)%id_atom(ii)%fpy(fp_index) = &
                        amat_features_ML(nn)%id_atom(ii)%fpy(fp_index) + tmpy
                    amat_features_ML(nn)%id_atom(ii)%fpz(fp_index) = &
                        amat_features_ML(nn)%id_atom(ii)%fpz(fp_index) + tmpz
                  end do ! three-body terms
                else if (i_species == k_species) then
                  ! X-AX
                  if (ia_glob==1 .and. flag_debug_mlff) then
                    write(*,104) param_start, param_end, &
                      i_species,j_species,k_species,amat(nn)%n_nab(ii), 'X-AX'
                  end if

                  do param_index= param_start, param_end
                    fp_index = param_index + shift_dim

                    ! Todo: check dimension of params_2b(j)
                    gx_index = param_index - param_start + 1
                    eta1 = descriptor_params(i_species)%params_3b(species_order3b)%eta1(gx_index)
                    eta2 = descriptor_params(i_species)%params_3b(species_order3b)%eta2(gx_index)
                    eta3 = descriptor_params(i_species)%params_3b(species_order3b)%eta3(gx_index)

                    ! transfter for efficiency
                    eta1 = 1.0/eta1 ** 2
                    eta2 = 1.0/eta2 ** 2
                    eta3 = 1.0/eta3 ** 2

                    ! Terms for A-A-A
                    tmp123 = exp(-eta1 * rik0 ** 2 - eta2 * rij0 ** 2 - eta3 * rjk0 ** 2)
                    tmp132 = exp(-eta1 * rik0 ** 2 - eta3 * rij0 ** 2 - eta2 * rjk0 ** 2)

                    proj_x = 0.0
                    proj_y = 0.0
                    proj_z = 0.0
                    proj_x = proj_x + tmp123 * (xik0 * eta1 + xij0 * eta2)
                    proj_x = proj_x + tmp132 * (xik0 * eta1 + xij0 * eta3)

                    proj_y = proj_y + tmp123 * (yik0 * eta1 + yij0 * eta2)
                    proj_y = proj_y + tmp132 * (yik0 * eta1 + yij0 * eta3)

                    proj_z = proj_z + tmp123 * (zik0 * eta1 + zij0 * eta2)
                    proj_z = proj_z + tmp132 * (zik0 * eta1 + zij0 * eta3)

                    tmpx = frc_3b * proj_x
                    tmpy = frc_3b * proj_y
                    tmpz = frc_3b * proj_z

                    amat_features_ML(nn)%id_atom(ii)%fpx(fp_index) = &
                        amat_features_ML(nn)%id_atom(ii)%fpx(fp_index) + tmpx
                    amat_features_ML(nn)%id_atom(ii)%fpy(fp_index) = &
                        amat_features_ML(nn)%id_atom(ii)%fpy(fp_index) + tmpy
                    amat_features_ML(nn)%id_atom(ii)%fpz(fp_index) = &
                        amat_features_ML(nn)%id_atom(ii)%fpz(fp_index) + tmpz
                  end do ! three-body terms
                else if (j_species < k_species) then
                  !! X-AB
                  if (ia_glob==1 .and. flag_debug_mlff) then
                  write(*,104) param_start, param_end, &
                      i_species,j_species,k_species,amat(nn)%n_nab(ii), 'X-AB'
                  end if

                  do param_index= param_start, param_end
                    fp_index = param_index + shift_dim

                    ! Todo: check dimension of params_2b(j)
                    gx_index = param_index - param_start + 1
                    eta1 = descriptor_params(i_species)%params_3b(species_order3b)%eta1(gx_index)
                    eta2 = descriptor_params(i_species)%params_3b(species_order3b)%eta2(gx_index)
                    eta3 = descriptor_params(i_species)%params_3b(species_order3b)%eta3(gx_index)

                    ! transfter for efficiency
                    eta1 = 1.0/eta1 ** 2
                    eta2 = 1.0/eta2 ** 2
                    eta3 = 1.0/eta3 ** 2

                    ! Terms for A-AB or B-AB
                    tmp123 = exp(-eta1 * rij0 ** 2 - eta2 * rik0 ** 2 - eta3 * rjk0 ** 2)

                    proj_x = 0.0
                    proj_y = 0.0
                    proj_z = 0.0
                    proj_x = proj_x + tmp123 * (xij0 * eta1 + xik0 * eta2)
                    proj_y = proj_y + tmp123 * (yij0 * eta1 + yik0 * eta2)
                    proj_z = proj_z + tmp123 * (zij0 * eta1 + zik0 * eta2)

                    tmpx = frc_3b * proj_x
                    tmpy = frc_3b * proj_y
                    tmpz = frc_3b * proj_z

                    amat_features_ML(nn)%id_atom(ii)%fpx(fp_index) = &
                        amat_features_ML(nn)%id_atom(ii)%fpx(fp_index) + tmpx
                    amat_features_ML(nn)%id_atom(ii)%fpy(fp_index) = &
                        amat_features_ML(nn)%id_atom(ii)%fpy(fp_index) + tmpy
                    amat_features_ML(nn)%id_atom(ii)%fpz(fp_index) = &
                        amat_features_ML(nn)%id_atom(ii)%fpz(fp_index) + tmpz
                  end do ! three-body terms
                else
                  !! X-BA
                  if (ia_glob==1 .and. flag_debug_mlff) then
                  write(*,104) param_start, param_end, &
                      i_species,j_species,k_species,amat(nn)%n_nab(ii), 'X-AB'
                  end if

                  do param_index= param_start, param_end
                    fp_index = param_index + shift_dim

                    ! Todo: check dimension of params_2b(j)
                    gx_index = param_index - param_start + 1
                    eta1 = descriptor_params(i_species)%params_3b(species_order3b)%eta1(gx_index)
                    eta2 = descriptor_params(i_species)%params_3b(species_order3b)%eta2(gx_index)
                    eta3 = descriptor_params(i_species)%params_3b(species_order3b)%eta3(gx_index)

                    ! transfter for efficiency
                    eta1 = 1.0/eta1 ** 2
                    eta2 = 1.0/eta2 ** 2
                    eta3 = 1.0/eta3 ** 2

                    ! Terms for A-AB or B-AB
                    tmp123 = exp(-eta1 * rik0 ** 2 - eta2 * rij0 ** 2 - eta3 * rjk0 ** 2)

                    proj_x = 0.0
                    proj_y = 0.0
                    proj_z = 0.0
                    proj_x = proj_x + tmp123 * (xik0 * eta1 + xij0 * eta2)
                    proj_y = proj_y + tmp123 * (yik0 * eta1 + yij0 * eta2)
                    proj_z = proj_z + tmp123 * (zik0 * eta1 + zij0 * eta2)

                    tmpx = frc_3b * proj_x
                    tmpy = frc_3b * proj_y
                    tmpz = frc_3b * proj_z

                    amat_features_ML(nn)%id_atom(ii)%fpx(fp_index) = &
                        amat_features_ML(nn)%id_atom(ii)%fpx(fp_index) + tmpx
                    amat_features_ML(nn)%id_atom(ii)%fpy(fp_index) = &
                        amat_features_ML(nn)%id_atom(ii)%fpy(fp_index) + tmpy
                    amat_features_ML(nn)%id_atom(ii)%fpz(fp_index) = &
                        amat_features_ML(nn)%id_atom(ii)%fpz(fp_index) + tmpz
                  end do ! three-body terms
                end if ! X-XA, X-AX, X-AB, X-BA
              end if

            end do ! k, three-body
            !write(*,*) 'after part_nd_nabs in get_feature_split', inode
          end do ! j, two-body
          inp=inp+1  ! Indexes primary-set atoms
        enddo ! End ni in prim%nm_nodgroup
      else
        ! if no atoms, do nothing, and turn to next partition
        if(flag_debug_mlff) &
            write(*, *) 'Warning: No atoms in this partition get_feature_split', inode, nn
      endif ! End if(prim%nm_nodgroup>0)
    enddo ! End part_on_node
    102 format('check fp 2b: start end ',2i5,' species= ',3i3, ' num_nb', i5)
    103 format('check fp 3b: ',5i5,' species= ',4i3, ' num_nb', i5)
    104 format('check fp 3b: start end ',2i5,' species= ',3i3, ' num_nb', i5,' ',a12)
    201 format('before',4e16.6)
    202 format('after',4e16.6)
  end subroutine get_feature_split
!!***

!!****f* mlff_type_module/clean_up_ml *
!!
!!  NAME
!!   clean_up_ml
!!  USAGE
!!
!!  PURPOSE
!!   Creates neighbour lists for a matrix range rcut
!!  INPUTS
!!   type(primary_set) :: prim
!!   type(cover_set) :: gcs
!!   type(matrix), dimension(:) :: amat
!!   real(double) :: rcut
!!  USES
!!
!!  AUTHOR
!!   Jianbo Lin
!!  CREATION DATE
!!   2022/07/20
!!  MODIFICATION HISTORY
!!   2022/07/27 J.Lin
!!    Modified to get pair information for machine learing
!!  SOURCE
!!
  subroutine clean_up_ml_acsf(amat,amat_features_ML,descriptor_params)
    implicit none

    ! Passed variables
    type(matrix_ML), dimension (:), allocatable, intent(out) :: amat
    type(features_ML), dimension(:), allocatable, intent(out) :: amat_features_ML
    type(acsf_param), dimension(:), allocatable, intent(out) :: descriptor_params

  end subroutine clean_up_ml_acsf
!!***

!!****f* mlff_type_module/clean_up_ml_split *
!!
!!  NAME
!!   clean_up_ml_split
!!  USAGE
!!
!!  PURPOSE
!!   Creates neighbour lists for a matrix range rcut
!!  INPUTS
!!   type(primary_set) :: prim
!!   type(cover_set) :: gcs
!!   type(matrix), dimension(:) :: amat
!!   real(double) :: rcut
!!  USES
!!
!!  AUTHOR
!!   Jianbo Lin
!!  CREATION DATE
!!   2022/07/20
!!  MODIFICATION HISTORY
!!   2022/07/27 J.Lin
!!    Modified to get pair information for machine learing
!!  SOURCE
!!
  subroutine clean_up_ml_split(amat,amat_features_ML,descriptor_params)
    implicit none

    ! Passed variables
    type(matrix_ML), dimension (:), allocatable, intent(out) :: amat
    type(features_ML), dimension(:), allocatable, intent(out) :: amat_features_ML
    type(split_param), dimension(:), allocatable, intent(out) :: descriptor_params

  end subroutine clean_up_ml_split
!!***

!!****f* mlff_type_module/double_swap *
!!
!!  NAME
!!   double_swap
!!  USAGE
!!
!!  PURPOSE
!!   Swap two real numbers
!!  INPUTS
!!
!!  USES
!!
!!  AUTHOR
!!   Jianbo Lin
!!  CREATION DATE
!!   2023/07/20
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine double_swap(a,b)
    implicit none
    real(double) :: a,b,tmp

    tmp=a
    a=b
    b=tmp
  end subroutine double_swap
!!***

!!****f* mlff_type_module/int_swap *
!!
!!  NAME
!!   int_swap
!!  USAGE
!!
!!  PURPOSE
!!   Swap two integer numbers
!!  INPUTS
!!
!!  USES
!!
!!  AUTHOR
!!   Jianbo Lin
!!  CREATION DATE
!!   2023/07/20
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine int_swap(a,b)
    implicit none
    integer :: a,b,tmp

    tmp=a
    a=b
    b=tmp
  end subroutine int_swap
!!***

!!****f* mlff_type_module/split_string *
!!
!!  NAME
!!   split_string
!!  USAGE
!!
!!  PURPOSE
!!   Split string depending on delimeter to an array of string
!!  Useful for collecting element names in one line
!!  INPUTS
!!
!!  USES
!!
!!  AUTHOR
!!   Jianbo Lin
!!  CREATION DATE
!!   2023/08/24
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine split_string(string_in, string_out, delim, n_species)
    use GenComms, ONLY: inode, ionode

    implicit none
    character(len=180)             :: string_in
    character(len=30)              :: delim
    character(len=2), dimension(:) :: string_out
    integer :: n_species

    integer :: index, ii
    character(len=180)             :: string_tmp

    string_tmp = TRIM(adjustl(string_in))
    index = SCAN(string_tmp,delim)
    ii=1
    do while ((index .gt. 0) .and. (LEN(TRIM(string_tmp)) .gt.2))
      string_out(ii) = string_tmp(1:index-1)
      string_tmp = string_tmp(index+1:)
      string_tmp = TRIM(adjustl(string_tmp))

      ii=ii+1
      index = SCAN(string_in,delim)
    end do
    string_out(ii)=TRIM(string_tmp)

  end subroutine split_string
end module mlff_type
