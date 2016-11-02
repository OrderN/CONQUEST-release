! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module Hilbert3D
! ------------------------------------------------------------------------------
! Code area 1: initialisation
! ------------------------------------------------------------------------------

!****h* Conquest/Hilbert3D
! NAME
!   Hilbert3D
! PURPOSE
!   Generates a Hilbert curve for a box of dimension 2**ax by 2**ay by
!   2**az, where ax, ay and az can be different non-negative integers
! AUTHOR
!   L.Tong
! CREATION DATE
!   2012/11/22
! MODIFICATION HISTORY
!   2016/10/06 dave
!    Changed allocatable types in unpacked to pointers for F90 compatibility
! SOURCE
!
module Hilbert3D

  implicit none
  
  save
  private
  
  public ::                   &
       Hilbert3D_Initialise,  &
       Hilbert3D_IntToCoords, &
       Hilbert3D_CoordsToInt

  character(len=80) :: &
       RCSid = "$Id$"

  integer, dimension(3) :: Alpha
  integer, dimension(3) :: Direction
  integer, dimension(3) :: XYZ
  integer :: DimLine, DimSquare, DimCube
  integer :: NSquaresInLine, NCubesInSquare, NPointsInCube

  ! unpacked index and coordinate structures
  type unpacked
     integer :: N_l_levels, N_p_levels, N_c_levels
     integer, dimension(:), pointer :: l
     integer, dimension(:), pointer :: p
     integer, dimension(:), pointer :: c
  end type unpacked

contains

  !****f* Hilbert3D/InitialiseHilbert
  ! PURPOSE
  !   Initialises Hilbert3D
  ! USAGE
  !   call InitialiseHilbert(ax, ay, az)
  ! INPUTS
  !   integer ax: x dimension will have 2**ax grid points
  !   integer ay: y dimension will have 2**ax grid points
  !   integer az: z dimension will have 2**az grid points
  ! OUTPUT
  !   Sets and initialises internal parameters used by Hilbert3D
  ! AUTHOR
  !   L.Tong
  ! CREATION DATE
  !   2012/11/22
  ! MODIFICATION HISTORY
  ! SOURCE
  !
  subroutine Hilbert3D_Initialise(ax, ay, az)
    implicit none
    ! Passed variables
    integer, intent(in) :: ax, ay, az
    ! Local variables
    integer :: n, newn, temp, ii
    Alpha = (/ax, ay, az/)
    Direction = (/4, 2, 1/)
    XYZ = (/1, 2, 3/)
    ! sort Alpha so that the longest dimension is left-most, and
    ! shortest right-most
    n = 3
    do while (n > 1)
       newn = 1
       do ii = 2, n
          if (Alpha(ii-1) < Alpha(ii)) then
             ! swap Alpha
             temp = Alpha(ii)
             Alpha(ii) = Alpha(ii-1)
             Alpha(ii-1) = temp
             ! swap Direction
             temp = Direction(ii)
             Direction(ii) = Direction(ii-1)
             Direction(ii-1) = temp
             ! swap XYZ
             temp = XYZ(ii)
             XYZ(ii) = XYZ(ii-1)
             XYZ(ii-1) = temp
             ! set newn
             newn = ii
          end if
       end do
       n = newn
    end do
    ! Set the dimensions
    DimLine = 2**(Alpha(1) - Alpha(2))
    DimSquare = 2**(Alpha(2) - Alpha(3))
    DimCube = 2**Alpha(3)
    NSquaresInLine = DimLine
    NCubesInSquare = DimSquare**2
    NPointsInCube = DimCube**3
  end subroutine Hilbert3D_Initialise
  !*****

  !****f* Hilbert3D/Hilbert3D_IntToCoords
  ! PURPOSE
  !   Given Hilbert curve step index (counting from 1) as input,
  !   outputs the cartesian grid coordinates of the corresponding
  !   point in the Hilbert curve.
  ! USAGE
  !   call Hilbert3D_IntToCoords(ind, coords)
  ! INPUTS
  !   integer ind: input index
  ! OUTPUT
  !   integer coords(3): output coordinates
  ! AUTHOR
  !   L.Tong
  ! CREATION DATE
  !   2012/11/22
  ! MODIFICATION HISTORY
  ! SOURCE
  !
  subroutine Hilbert3D_IntToCoords(ind, coords)
    implicit none
    ! Passed parameters
    integer, intent(in) :: ind
    integer, dimension(:), intent(out) :: coords
    ! Local variables
    type(unpacked) :: unpkd_ind, unpkd_coords
    integer :: index, c_start, c_end, start, end, ii
    integer :: l_ind, p_ind, c_ind
    ! ind goes starts from 1, transform to index which starts from 0
    index = max(0, ind - 1)
    ! unpack ind, unpkd_ind is allocated by UnpackIndex
    call UnpackIndex(index, unpkd_ind)
    ! allocate the unpacked coordinates with same dimension as unpkd_ind
    call AllocUnpacked(unpkd_coords, unpkd_ind%N_l_levels, &
                       unpkd_ind%N_p_levels, unpkd_ind%N_c_levels)
    ! zero the l, p and c levels for unpkd_coords
    unpkd_coords%l = 0
    unpkd_coords%p = 0
    unpkd_coords%c = 0
    ! move along the outer-mose level of the line
    start = 0
    end = 1
    if (DimLine > 1) then
       ! recursively generate the 1D fractal coordinates
       do ii = 1, unpkd_coords%N_l_levels
          l_ind = unpkd_ind%l(ii)
          unpkd_coords%l(ii) = GrayEncodeTravel(start, end, 1, l_ind)
          call ChildStartEnd(start, end, 1, l_ind, c_start, c_end)
          start = c_start
          end = c_end
       end do
    end if
    ! move in the middle level squares, the highest level unit square
    ! always move along the direction of the higher level 1D line
    start = start * 2 ! need to convert from 1 to 10 (1D to 2D)
    end = end * 2
    if (DimSquare > 1) then
       do ii = 1, unpkd_coords%N_p_levels
          p_ind = unpkd_ind%p(ii)
          unpkd_coords%p(ii) = GrayEncodeTravel(start, end, 2, p_ind)
          call ChildStartEnd(start, end, 2, p_ind, c_start, c_end)
          start = c_start
          end = c_end
       end do
    end if
    ! move in the inner level cubes, the highest level unit cube
    ! follows the direction of the children of the upper level unit
    ! square, unless dimension of 2D curve is 1, in which case it
    ! follows the direction of the line
    start = start * 2 ! need to convert from 10 to 100 (2D to 3D)
    end = end * 2
    if (DimCube > 1) then
       do ii = 1, unpkd_coords%N_c_levels
          c_ind = unpkd_ind%c(ii)
          unpkd_coords%c(ii) = GrayEncodeTravel(start, end, 3, c_ind)
          call ChildStartEnd(start, end, 3, c_ind, c_start, c_end)
          start = c_start
          end = c_end
       end do
    end if
    ! pack unpkd_coords to give the cartesian coordinates
    call PackCoords(unpkd_coords, coords)
    ! deallocate memory
    call DeallocUnpacked(unpkd_ind)
    call DeallocUnpacked(unpkd_coords)
  end subroutine Hilbert3D_IntToCoords
  !*****
  
  !****f* Hilbert3D/Hilbert3D_CoordsToInt
  ! PURPOSE
  !   Inverse transorm of Hilbert3D_IntToCoords.
  !   
  !   Given a cartesian grid coordinate (nx,ny,nz), returns the
  !   corresponding step index (counting from 1) of the Hilbert Curve.
  ! USAGE
  !   call Hilbert3D_CoordsToInt(coords, ind)
  ! INPUTS
  !   integer coords: input cartesian coordinate
  ! OUTPUT
  !   integer ind: step index (counting from 1) of coords in Hilbert
  !                curve
  ! AUTHOR
  !   L.Tong
  ! CREATION DATE
  !   2012/11/25
  ! MODIFICATION HISTORY
  ! SOURCE
  !
  subroutine Hilbert3D_CoordsToInt(coords, ind)
    implicit none
    ! Passed parameters
    integer, dimension(:), intent(in) :: coords
    integer, intent(out) :: ind
    ! Local variables
    type(unpacked) :: unpkd_coords, unpkd_ind
    integer :: index, start, end, c_start, c_end, ii
    integer :: l_coords, p_coords, c_coords
    call UnpackCoords(coords, unpkd_coords)
    call AllocUnpacked(unpkd_ind, unpkd_coords%N_l_levels, &
                       unpkd_coords%N_p_levels, unpkd_coords%N_c_levels)
    ! zero the l, p and c levels for unpkd_ind
    unpkd_ind%l = 0
    unpkd_ind%p = 0
    unpkd_ind%c = 0
    ! move along the outer level 1D Hilbert curve
    start = 0
    end = 1
    if (DimLine > 1) then
       do ii = 1, unpkd_coords%N_l_levels
          l_coords = unpkd_coords%l(ii)
          unpkd_ind%l(ii) = GrayDecodeTravel(start, end, 1, l_coords)
          call ChildStartEnd(start, end, 1, unpkd_ind%l(ii), c_start, c_end)
          start = c_start
          end = c_end
       end do
    end if
    ! move in the middle level 2D Hilbert curve, the highest level
    ! square always has travel along direction of the line
    start = start * 2 ! need to convert from 1 to 10 (1D to 2D)
    end = end * 2
    if (DimSquare > 1) then
       do ii = 1, unpkd_coords%N_p_levels
          p_coords = unpkd_coords%p(ii)
          unpkd_ind%p(ii) = GrayDecodeTravel(start, end, 2, p_coords)
          call ChildStartEnd(start, end, 2, unpkd_ind%p(ii), c_start, c_end)
          start = c_start
          end = c_end
       end do
    end if
    ! move in the lower level 3D Hilbert curve, the highest level cube
    ! follows the direction of children of the 2D curve, unless
    ! dimension of 2D curve is 1, in which case it follows the
    ! direction of the 1D curve
    start = start * 2 ! need to convert from 10 to 100 (2D to 3D)
    end = end * 2
    if (DimCube > 1) then
       do ii = 1, unpkd_coords%N_c_levels
          c_coords = unpkd_coords%c(ii)
          unpkd_ind%c(ii) = GrayDecodeTravel(start, end, 3, c_coords)
          call ChildStartEnd(start, end, 3, unpkd_ind%c(ii), c_start, c_end)
          start = c_start
          end = c_end
       end do
    end if
    call PackIndex(unpkd_ind, index)
    ! index starts counting from 0, so correct for ind which counts
    ! from 1
    ind = index + 1
    ! deallocate memory
    call DeallocUnpacked(unpkd_ind)
    call DeallocUnpacked(unpkd_coords)
  end subroutine Hilbert3D_CoordsToInt
  !*****
  
  !****f* Hilbert3D/AllocUnpack
  ! PURPOSE
  !   Allocates the unpacked index or coordinates
  ! USAGE
  !   call AllocUnpack(unpkd, nLLevels, nPLevels, nClevels)
  ! INPUTS
  !   integer nLLevels: number of fractal levels in the 1D Hilbert curve
  !   integer nPLevels: number of fractal levels in the 2D Hilbert curve
  !   integer nCLevels: number of fractal levels in the 3D Hilbert curve
  ! OUTPUT
  !   type(unpacked) unpkd: the unpacked type to be allocated
  ! RETURN VALUE
  ! AUTHOR
  !   L.Tong
  ! CREATION DATE
  !   2012/11/22
  ! MODIFICATION HISTORY
  ! SOURCE
  !
  subroutine AllocUnpacked(unpkd, nLLevels, nPLevels, nCLevels)
    use global_module, only: area_init
    use GenComms,      only: cq_abort
    use memory_module, only: reg_alloc_mem, type_int
    implicit none
    ! Passed parameters
    integer, intent(in) :: nLLevels, nPLevels, nCLevels
    type(unpacked), intent(out) :: unpkd
    ! Local variable
    integer :: stat
    allocate(unpkd%l(nLLevels), unpkd%p(nPLevels), unpkd%c(nCLevels), &
             STAT=stat)
    if (stat /= 0) then
         call cq_abort("AllocUnpack: Error allocating unpkd. ", stat) 
    end if
    call reg_alloc_mem(area_init, nLLevels + nPLevels + nCLevels, &
                       type_int)
    unpkd%N_l_levels = nLLevels
    unpkd%N_p_levels = nPLevels
    unpkd%N_c_levels = nCLevels
  end subroutine AllocUnpacked
  !*****

  !****f* Hilbert3D/DeallocUnpack(unpkd)
  ! PURPOSE
  !   Deallocates a unpacked type
  ! USAGE
  !   call DeallocUnpack(unpkd)
  ! INPUTS
  !   type(unpacked) unpkd: unpacked type to be deallocated
  ! OUTPUT
  !   unpkd deallocated
  ! AUTHOR
  !   L.Tong
  ! CREATION DATE
  !   2012/11/22
  ! MODIFICATION HISTORY
  ! SOURCE
  !
  subroutine DeallocUnpacked(unpkd)
    use global_module, only: area_init
    use GenComms,      only: cq_abort
    use memory_module, only: reg_dealloc_mem, type_int
    implicit none
    ! Passed parameters
    type(unpacked), intent(inout) :: unpkd
    ! Local variable
    integer :: stat
    deallocate(unpkd%l, unpkd%p, unpkd%c, STAT=stat)
    if (stat /= 0) then
       call cq_abort("DeallocUnpack: Error deallocating unpkd. ", stat)
    end if
    call reg_dealloc_mem(area_init, &
                         unpkd%N_l_levels + &
                         unpkd%N_p_levels + &
                         unpkd%N_c_levels, &
                         type_int)
    unpkd%N_l_levels = 0
    unpkd%N_p_levels = 0
    unpkd%N_c_levels = 0
  end subroutine DeallocUnpacked
  !*****

  !****f* Hilbert3D/UnpackIndex
  ! PURPOSE
  !   Unpacks the step index (start counting from 0) in Hilbert curve
  !   into the decomposed set of indices in the corresponding
  !   different levels of fractal unit cubes and the 1D, 2D and 3D
  !   Hilbert curves
  ! USAGE
  !   call UnpackIndex(ind, unpkdInd)
  ! INPUTS
  !   integer ind: step index which starts counting from 0
  ! OUTPUT
  !   integer unpkdInd: unpacked index (allocated in the process)
  ! AUTHOR
  !   L.Tong
  ! CREATION DATE
  !   2012/11/22
  ! MODIFICATION HISTORY
  ! SOURCE
  !
  subroutine UnpackIndex(ind, unpkdInd)
    implicit none
    ! Passed parameters
    integer, intent(in) :: ind
    type(unpacked), intent(out) :: unpkdInd
    ! Local variables
    integer :: ic, ip, l, p, c, ii
    integer :: n_l_levels, n_p_levels, n_c_levels
    ic = ind / NPointsInCube
    c = mod(ind, NPointsInCube)
    ip = ic / NCubesInSquare
    p = mod(ic, NCubesInSquare)
    l = ip
    ! work out the number of fractal levels in the line, square and cube
    n_l_levels = MinLen(NSquaresInLine-1, 2)  ! 2 points in unit line
    n_p_levels = MinLen(NCubesInSquare-1, 4)  ! 4 points in unit square
    n_c_levels = MinLen(NPointsInCube-1, 8)   ! 8 points in unit cube
    ! allocate unpkd_ind
    call AllocUnpacked(unpkdInd, n_l_levels, n_p_levels, n_c_levels)
    ! unpack the indices corresponding to l of outer line
    do ii = n_l_levels, 1, -1
       unpkdInd%l(ii) = mod(l, 2)
       l = l / 2
    end do
    ! unpack the indices corresponding to p of middle squares
    do ii = n_p_levels, 1, -1
       unpkdInd%p(ii) = mod(p, 4)
       p = p / 4
    end do
    ! unpack the indices corresponding to c of inner cubes
    do ii = n_c_levels, 1, -1
       unpkdInd%c(ii) = mod(c, 8)
       c = c / 8
    end do
  end subroutine UnpackIndex
  !*****
  
  !****f* Hilbert3D/PackIndex
  ! PURPOSE
  !   Packs a type(unpacked) index into an integer step index
  !   (counting starts from 0) of the Hilbert curve. This is an inverse
  !   transformation of UnpackIndex
  ! USAGE
  !   call PackIndex(unpkdInd, ind)
  ! INPUTS
  !   type(unpacked) unpkdInd: unpacked index
  ! OUTPUT
  !   integer ind:  packed index
  ! AUTHOR
  !   L.Tong
  ! CREATION DATE
  !   2012/11/22
  ! MODIFICATION HISTORY
  ! SOURCE
  !
  subroutine PackIndex(unpkdInd, ind)
    implicit none
    ! Passed parameters
    type(unpacked), intent(in) :: unpkdInd
    integer, intent(out) :: ind
    ! Local variables
    integer :: l, p, c, ip, ic, ii
    ! pack the indices corresponding to outer lines
    l = unpkdInd%l(1)
    do ii = 2, unpkdInd%N_l_levels
       l = l * 2 + unpkdInd%l(ii)
    end do
    ! pack the indices corresponding to middle squares
    p = unpkdInd%p(1)
    do ii = 2, unpkdInd%N_p_levels
       p = p * 4 + unpkdInd%p(ii)
    end do
    ! pack the indices corresponding to inndr cubes
    c = unpkdInd%c(1)
    do ii = 2, unpkdInd%N_c_levels
       c = c * 8 + unpkdInd%c(ii)
    end do
    ! pack l, p and c into a single index
    ip = l
    ic = ip * NCubesInSquare + p
    ind = ic * NPointsInCube + c
  end subroutine PackIndex
  !*****
  
  !****f* Hilbert3D/UnpackCoords
  ! PURPOSE
  !   Given cartesian coordinates, decomposes the coordinate into
  !   the corresponding binary code coordinates of the different
  !   levels of unit cubes and the 1D, 2D and 3D Hilbert curves.
  !
  !   Note the binary coordinates for unit cube are:
  !     for 1D is 0, 1
  !     for 2D is 00, 10, 01, 11
  !     for 3D is 000, 100, 010, 001, 110, 101, 011, 111
  !   And the left most bit always corresponds to the longgest
  !   dimension and the right most the shortest dimension. The
  !   required rotation from the input (x,y,z) in cartesian format
  !   are automatically taken care of.
  ! USAGE
  !   call UnpackCoords(coords, unpkdCoords)
  ! INPUTS
  !   integer coords(3): input cartesian coordinates
  ! OUTPUT
  !   type(unpacked) unpkdCoords: output unpacked coordinates
  ! AUTHOR
  !   L.Tong
  ! CREATION DATE
  !   2012/11/22
  ! MODIFICATION HISTORY
  ! SOURCE
  !
  subroutine UnpackCoords(coords, unpkdCoords)
    implicit none
    ! Passed parameter
    integer, dimension(3), intent(in) :: coords
    type(unpacked), intent(out) :: unpkdCoords
    ! Local variable
    integer :: ip, ii
    integer :: n_l_levels, n_p_levels, n_c_levels
    integer, dimension(1) :: l
    integer, dimension(2) :: p, ic
    integer, dimension(3) :: c
    do ii = 1, 2
       ic(ii) = coords(XYZ(ii)) / DimCube
    end do
    do ii = 1, 3
       c(ii) = mod(coords(XYZ(ii)), DimCube)
    end do
    ip = ic(1) / DimSquare
    do ii = 1, 2
       p(ii) = mod(ic(ii), DimSquare)
    end do
    l(1) = ip
    ! allocate unpkdCoords
    n_l_levels = MinLen(DimLine-1, 2)
    n_p_levels = MinLen(DimSquare-1, 2)
    n_c_levels = MinLen(DimCube-1, 2)
    call AllocUnpacked(unpkdCoords, n_l_levels, n_p_levels, &
                       n_c_levels)
    ! unpack coordinates corresponding to l, p, c
    call TransposeBits(1, l, n_l_levels, unpkdCoords%l)
    call TransposeBits(2, p, n_p_levels, unpkdCoords%p)
    call TransposeBits(3, c, n_c_levels, unpkdCoords%c)
  end subroutine UnpackCoords
  !*****
  
  !****f* Hilbert3D/PackCoords
  ! PURPOSE
  !   Inverse transformation of UnpackCoords
  ! USAGE
  !   call PackCoords(unpkdCoords, coords)
  ! INPUTS
  !   type(unpacked) unpkdCoords: input unpacked coordinates
  ! OUTPUT
  !   integer coords(3): output cartesian coordinate
  ! AUTHOR
  !   L.Tong
  ! CREATION DATE
  !   2012/11/22
  ! MODIFICATION HISTORY
  ! SOURCE
  !
  subroutine PackCoords(unpkdCoords, coords)
    implicit none
    ! Passed parameters
    type(unpacked), intent(in) :: unpkdCoords
    integer, dimension(:), intent(out) :: coords
    ! Local variables
    integer, dimension(1) :: l
    integer, dimension(2) :: p, ic
    integer, dimension(3) :: c
    integer :: ip
    call TransposeBits(unpkdCoords%N_l_levels, unpkdCoords%l, 1, l)
    call TransposeBits(unpkdCoords%N_p_levels, unpkdCoords%p, 2, p)
    call TransposeBits(unpkdCoords%N_c_levels, unpkdCoords%c, 3, c)
    ip = l(1)
    ic(1) = ip * DimSquare + p(1)
    ic(2) = p(2)
    coords(XYZ(1)) = ic(1) * DimCube + c(1)
    coords(XYZ(2)) = ic(2) * DimCube + c(2)
    coords(XYZ(3)) = c(3)
  end subroutine PackCoords
  !*****

  !****f* Hilbert3D/ChildStartEnd
  ! PURPOSE
  !   Given parent starting position and ending position of a Gray
  !   code sequence, and the parent step index, works out child 
  !   starting position and ending position for the fractal expansion
  !   of a canonical Hilbert curve.
  ! USAGE
  !   call ChildStartEnd(parentStart, parentEnd, nDims, parentInd,
  !                      childStart, childEnd)
  ! INPUTS
  !   integer parentStart: starting Gray code in parent sequence
  !   integer parentEnd:   ending Gray code in parent sequence
  !   integer nDims:       number of dimensions for Hilbert curve,
  !                        equal to the number of bits for Gray code
  !   integer parentInd:   step index in parent sequence, index
  !                        starts counting from 0
  ! OUTPUT
  !   integer childStart:  starting Gray code in child sequence
  !   integer childEnd:    ending Gray code in child sequence
  ! AUTHOR
  !   L.Tong
  ! CREATION DATE
  !   2012/11/25
  ! MODIFICATION HISTORY
  ! SOURCE
  !
  subroutine ChildStartEnd(parentStart, parentEnd, nDims, parentInd, &
                           childStart, childEnd)
    implicit none
    ! Passed parameters
    integer, intent(in) :: parentStart, parentEnd
    integer, intent(in) :: nDims, parentInd
    integer, intent(out) :: childStart, childEnd
    ! Local variables
    integer :: maximum, start_ind, end_ind
    ! maximum value or the nDim number of bits
    maximum = 2**nDims - 1
    ! next lower even number or 0
    start_ind = max(0, iand((parentInd - 1), not(1)))
    ! next higher odd number or maximum
    end_ind = min(maximum, ior((parentInd + 1), 1))
    childStart = GrayEncodeTravel(parentStart, parentEnd, nDims, start_ind)
    childEnd = GrayEncodeTravel(parentStart, parentEnd, nDims, end_ind)
  end subroutine ChildStartEnd
  !*****
  
  !****f* Hilbert3D/GrayEncode
  ! PURPOSE
  !   Gray code encoder for canonical reflective Gray code sequence.
  !   Remember the sequence starts counting from 0
  ! USAGE
  !   ith_gray_code = GrayEncode(i)
  ! INPUTS
  !   integer i: input integer, read as a binary number, corresponds
  !              to i-th code in canonical Gray sequence
  ! RETURN VALUE
  !   integer GrayEncode: the gray code represented as an integer
  ! AUTHOR
  !   L.Tong
  ! CREATION DATE
  !   2012/11/22
  ! MODIFICATION HISTORY
  ! SOURCE
  !
  function GrayEncode(i)
    implicit none
    ! Passed variables
    integer, intent(in) :: i
    ! Return value
    integer :: GrayEncode
    GrayEncode = ieor(i, ishft(i,-1))
  end function GrayEncode
  !*****
  
  !****f* Hilbert3D/GrayDecode
  ! PURPOSE
  !   Inverse operation of GrayEncode. Given a canonical Gray code,
  !   returns the integer index of the code in the Gray code sequence.
  !   Remember the sequence starts counting from 0
  ! USAGE
  !   i = GrayEncode(gn)
  ! INPUTS
  !   integer gn: input Gray code in integer representation
  ! RETURN VALUE
  !   integer GrayDecode: the integer index of gn in Gray code sequence
  ! AUTHOR
  !   L.Tong
  ! CREATION DATE
  !   2012/11/25
  ! MODIFICATION HISTORY
  ! SOURCE
  !
  function GrayDecode(gn)
    implicit none
    ! Passed variables
    integer, intent(in) :: gn
    ! Return value
    integer :: GrayDecode
    ! Local variables
    integer :: bn, ii
    bn = gn
    ii = 1
    do while (ishft(gn, -ii) > 0)
       bn = ieor(bn, ishft(gn, -ii))
       ii = ii + 1
    end do
    GrayDecode = bn
  end function GrayDecode
  !*****
  
  !****f* Hilbert3D/GrayEncodeTravel
  ! PURPOSE
  !   Let n = nBits, the canonical n-bit Gray code starts from 0,
  !   and travels to 2**(n-1). This corresponds to the Hilbert curve
  !   travelling in the n-th bit direction (counting from the Right
  !   starting from 1) in a n-dimensional unit cube. Note that j-th
  !   bit corresponds to the j-th coordinate of the n-dimensional
  !   unit cube.
  !
  !   This function rotates the Hilbert curve in the unit cube of
  !   nBits dimension, so that the path starts from start and ends
  !   at end, and for a given integer i, function returns the
  !   corresponding transformed Gray code sequence element.
  !
  !   Note that the transformed Gray code sequence may no longer be
  !   canonical. Indeed the canonical reflective Gray code while
  !   cyclic, does not give all Hilbert curves with desired staring
  !   and ending points in the given unit cube. Also note that
  !   simply defining the starting and ending point in a n-dimension
  !   unit cube (n>2) does NOT uniqually define the Hilbert curve
  !   (or the corresponding Gray code sequence): For example in 3D:
  !
  !   000 001 011 010 110 111 101 100
  !   000 010 011 001 101 111 110 100
  !
  !   both are Gray code sequences that start from 000 and end at
  !   100. The first one is the canonical Gray code sequence, the
  !   second one is one that is conjugate to the first sequence.
  !
  !   This function will only produce ONE particular solution given
  !   the starting and ending positions, and the solution will NOT
  !   necessarily be a cyclic permutation of the canonical Gray code
  !   sequence, even if a cyclic permutation of the canonical
  !   sequence could give a Hilbert curve that satisfy sthe start
  !   and end conditions.
  !   
  !   It is assumed that the start and end must be such that the
  !   travel direction (start ^ end) must be only contain one
  !   none-zero bit. That is the curve must be a Hilbert curve, with
  !   overall travel in x, y, or z, or ... (>3D) directions.
  !
  !   Since the canonical Gray code travel is always 2**(n-1), to
  !   transform to a code sequence with travel equal to (start ^
  !   end) = 2**t, we must rotate the 'axis' i.e. the bits to the
  !   left by t+1 bits.
  !
  !   To rotate the bit ordering of an n-bit number by m bits to the
  !   left, we just need to call the fortran intrinsic ishftc:
  !
  !     rotated_b = ishftc(b, m, n-bit)
  !   
  !   However for the purpose of this function it may not be as
  !   efficient as doing
  !
  !   1. shift b to the left by m bits
  !   2. shift b to the right by n-m bits
  !   3. merge the 1. and 2. using, say | (bitwise OR) operator
  !   4. mask out the hang-over bits to the left by apply a mask,
  !      where mask = 1<<n - 1 = 2**n - 1 (i.e. n bits all equal to
  !      1). Then mask with & (bitwise AND) operation.
  !
  !   This is because in order to work out m = t+1, we need to take
  !   log of travel = start ^ end. But by doing the above 4 steps, we
  !   can shift b to the left simply by multiplying travel * 2 =
  !   travel << 1.
  !
  !   Once we obtained a Gray code sequence that starts from 0 and
  !   ends at (start ^ end), to transform to start point start, all
  !   we need to do is to ^ with start, as start ^ 0 = start, and
  !   start ^ start = 0. Note that bn ^ Gn (for any n-bit binary,
  !   and Gn any n-bit Gray code sequence) is still a Gray code
  !   sequence. This is because considering adjacent Gray codes in
  !   the sequence, only one bit is changed. For all the bits that
  !   has not been changed, XOR with an n-bit binary gives same
  !   results; for the bit that is fliped, the result of XOR
  !   operation will either flip or not-flip the bit, in either case
  !   they will be different.
  !
  !   Remember the sequence starts counting from 0
  ! USAGE
  !   gray_code = GrayEncodeTravel(start, end, nBits, ind)
  ! INPUTS
  !   integer start: starting gray code in the sequence in integer representation
  !   integer end:   ending gray code in the sequence in integer representation
  !   integer nBits: number of bits in each code in the sequence
  !   integer ind:   integer index of the particular code in the sequence
  ! RETURN VALUE
  !   integer GrayEndodeTravel:  The particular Gray code
  !                              corresponding to ind in the sequence 
  ! AUTHOR
  !   L.Tong
  ! CREATION DATE
  !   2012/11/22
  ! MODIFICATION HISTORY
  ! SOURCE
  !
  function GrayEncodeTravel(start, end, nBits, ind)
    implicit none
    ! Passed variables
    integer, intent(in) :: start, end, nBits, ind
    ! Returned variable
    integer :: GrayEncodeTravel
    ! Local variables
    integer :: travel, mask, gn, rotated_gn
    travel = ieor(start, end)
    mask = ishft(1, nBits) - 1
    gn = GrayEncode(ind) * ishft(travel, 1)
    rotated_gn = iand(ior(gn, ishft(gn, -nBits)), mask)
    GrayEncodeTravel = ieor(rotated_gn, start) 
  end function GrayEncodeTravel
  !*****
  
  !****f* Hilbert3D/GrayDecodeTravel
  ! PURPOSE
  !   Inverse transform of GrayEncode Travel.
  !   Remember the sequence starts counting from 0
  ! USAGE
  !   ind = GrayDecodeTravel(start, end, nBits, gn)
  ! INPUTS
  !   integer start: starting Gray code of the sequence
  !   integer end:   ending Gray code of the sequence
  !   integer nBits: number of bits in the Gray codes
  !   integer gn:    the particular gray code whose index in sequence
  !                  we wish to find
  ! RETURN VALUE
  !   integer GrayDecodeTravel: the index in sequence of gn
  ! AUTHOR
  !   L.Tong
  ! CREATION DATE
  !   2012/11/22
  ! MODIFICATION HISTORY
  ! SOURCE
  !
  function GrayDecodeTravel(start, end, nBits, gn)
    implicit none
    ! Passed variables
    integer, intent(in) :: start, end, nBits, gn
    ! Returned variable
    integer :: GrayDecodeTravel
    ! Local variable
    integer :: travel, mask, rg, rotated_gn
    travel = ieor(start, end)
    mask = ishft(1, nBits) - 1
    rg = ishft(ieor(gn, start), nBits) / ishft(travel, 1)
    rotated_gn = iand(ior(rg, ishft(rg, -nBits)), mask)
    GrayDecodeTravel = GrayDecode(rotated_gn)
  end function GrayDecodeTravel
  !*****

  !****f* Hilbert3D/MinLen
  ! PURPOSE
  !   Give the minimum number of digits required to provide a
  !   representation of a given integer in a given base
  ! USAGE
  !   len = MinLen(n, base)
  ! INPUTS
  !   integer n: input integer
  !   integer base: base or representation
  ! RETURN VALUE
  !   integer MinLem: minimum number of digits needed to represent n
  ! AUTHOR
  !   L.Tong
  ! CREATION DATE
  !   2012/11/25
  ! MODIFICATION HISTORY
  ! SOURCE
  !
  function MinLen(n, base)
    implicit none
    ! Passed parameters
    integer, intent(in) :: n, base
    ! Returned value
    integer :: MinLen
    ! local variable
    integer :: num
    num = n
    MinLen = 1
    do while (num >= base)
       num = num / base
       MinLen = MinLen + 1
    end do
  end function MinLen
  !*****
  
  !****f* Hilbert3D/TransposeBits(nSrc, src, nDest, dest)
  ! PURPOSE
  !   Given array src of binary numbers each nDest bits long,
  !   calculates an array of nDest binary numbers each with nSrc bits
  !   in dest, such that the values are results of matrix transpose
  !   when src is vieweed as a matrix of bits.
  !
  !   For example: src = (010, 110, 001, 100), nSrc = 4, nDest = 3
  !   then dest is (0101, 1100, 0010)
  !
  !   0 1 0       \    0 1 0 1
  !   1 1 0   -----\   1 1 0 0
  !   0 0 1   -----/   0 0 1 0
  !   1 0 0       /
  ! USAGE
  !   call TransposeBits(nSrc, src, nDest, dest)
  ! INPUTS
  !   integer nSrc:        number of binary numbers in src
  !   integer src(nSrc):   input array of binary numbers
  !   integer nDest:       number of bits in binary number in src
  ! OUTPUT
  !   integer dest(nDest): the transposed array
  ! AUTHOR
  !   L.Tong
  ! CREATION DATE
  !   2012/11/25
  ! MODIFICATION HISTORY
  ! SOURCE
  !
  subroutine TransposeBits(nSrc, src, nDest, dest)
    implicit none
    ! Passed parameters
    integer, intent(in) :: nSrc, nDest
    integer, dimension(:), intent(in) :: src
    integer, dimension(:), intent(out) :: dest
    ! Local variables
    ! src is going to be a small array, so just use automatic arrays
    integer, dimension(nSrc) :: input ! copy of src
    integer :: ii, jj, output
    input = src
    do ii = nDest, 1, -1
       output = 0
       do jj = 1, nSrc
          output = ishft(output, 1) + mod(input(jj), 2)
          input(jj) = ishft(input(jj), -1)
       end do
       dest(ii) = output
    end do
  end subroutine TransposeBits
  !*****
  
end module Hilbert3D
!*****
