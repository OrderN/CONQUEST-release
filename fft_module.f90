! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module fft_module
! ------------------------------------------------------------------------------
! Code area 5: Self-consistency and charge density
! ------------------------------------------------------------------------------

!!****h* Conquest/fft_module
!!  NAME
!!   fft_module
!!  PURPOSE
!!   Collects all variables and subroutines associated with
!!   ffts into one place
!!  USES
!!   common, datatypes, dimens, GenComms, grid_index, numbers
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   31/05/2001 dave
!!  MODIFICATION HISTORY
!!   01/06/2001 dave
!!    Added RCS Log tag
!!   04/06/2001 dave
!!    Removed the dependency on nodes = 2^n in set_fft_map (I think)
!!    by replacing all IEOR loops with MPI_alltoall calls.
!!   06/06/2001 dave
!!    Added target attribute to pr?? and pack?? and corrected problems
!!    in fft3 and set_fft_map
!!   07/06/2001 dave
!!    Changed to get nodes from GenComms
!!   07/06/2001 dave
!!    Updated set_fft_map to use exch instead of alltoall and 
!!    cq_abort instead of stop
!!   18/03/2002 dave
!!    Added a little to the header and a static tag for object file id
!!   19/06/2002 dave
!!    Changed to use GPFA
!!   14:08, 29/07/2003 drb 
!!    Rewrote dependences on n_nodes in common and nodes in GenComms to be 
!!    mx_node and numprocs from global_module - much more transparent.
!!   08:33, 2003/09/05 dave
!!    Fixed MAJOR bug (by me - whoops !) where the SAME set of trigonometric factors was used for x, y AND z
!!    Also commented out MAX_FFT_SIZE variable which (it turns out) isn't being used !
!!   11:02, 14/11/2005 drb 
!!    Added array (recip_vector) to store reciprocal space vectors for grid 
!!   15:54, 27/04/2007 drb 
!!    Changed recip_vector to be (n,3) for speed
!!   2008/02/04 17:30 dave
!!    Changed for output to file not stdout
!!  SOURCE
!!
module fft_module

  use datatypes
  use global_module, ONLY: io_lun

  implicit none

  save

  ! RCS tag for object file identification   
  character(len=80), private :: RCSid = "$Id$"

  ! Previously ffttable
!  integer, parameter :: MAX_FFT_SIZE=8*128
  
!  real(double), dimension(100+MAX_FFT_SIZE*8) :: tablex,tabley,tablez
!  complex(double_cplx), dimension(3*(100+MAX_FFT_SIZE*8)) :: tablex,tabley,tablez

  ! Previously map
  integer, allocatable, dimension(:), target :: prbx
  integer, allocatable, dimension(:), target :: prxb
  integer, allocatable, dimension(:), target :: prxy
  integer, allocatable, dimension(:), target :: pryx
  integer, allocatable, dimension(:), target :: pryz
  integer, allocatable, dimension(:), target :: przy
  integer, allocatable, dimension(:) :: rprbx
  integer, allocatable, dimension(:) :: rprxb
  integer, allocatable, dimension(:) :: rprxy
  integer, allocatable, dimension(:) :: rpryx
  integer, allocatable, dimension(:) :: rpryz
  integer, allocatable, dimension(:) :: rprzy
  integer, allocatable, dimension(:), target :: packbx
  integer, allocatable, dimension(:), target :: packxb
  integer, allocatable, dimension(:), target :: packxy
  integer, allocatable, dimension(:), target :: packyx
  integer, allocatable, dimension(:), target :: packyz
  integer, allocatable, dimension(:), target :: packzy
  integer, dimension(:), pointer :: unpackbx, unpackxb
  integer, dimension(:), pointer :: psendxb, psendbx
  integer, dimension(:), pointer :: unpackyx, unpackxy
  integer, dimension(:), pointer :: psendyx, psendxy
  integer, dimension(:), pointer :: unpackyz, unpackzy
  integer, dimension(:), pointer :: psendyz, psendzy
  integer, allocatable, dimension(:) :: x_columns_node, y_columns_node, &
       z_columns_node

  real(double), allocatable, dimension(:) :: hartree_factor
  real(double), allocatable, dimension(:,:) :: recip_vector

  real(double), dimension(:), allocatable :: trigs_x, trigs_y, trigs_z

  integer :: i0 ! Gamma point location
  integer, allocatable, dimension(:) :: kerker_list
  integer :: size_kl
!!***

contains

! -----------------------------------------------------------
! Subroutine fft3
! -----------------------------------------------------------

!!****f* fft_module/fft3 *
!!
!!  NAME 
!!   fft3
!!  USAGE
!!   fft3(input data, output complex FFT of data, direction)
!!   fft3(data, cdata, isign)
!!  PURPOSE
!!   Perform an FFT on the charge density data. Note that the
!!   forward transform (isign = -1) starts from domain data
!!   structure and ends in columns parallel to the z axis, and the
!!   reverse transform does the opposite - in other words, the data
!!   is NOT stored in the same format when in reciprocal space.
!!
!!   Note on allocating trigs_x, y and z: 
!!   Note that the size of trigs_? is at least 2*(2^p+3^q+5^r), where n_grid_? factors 
!!   as 2^p*3^q*5^r.  If we set the size of trigs_? to 2*n_grid_? we'll be safe.
!!  INPUTS
!!   real(double), dimension(N_GRID_MAX) :: data
!!   complex(double_cplx), dimension(N_GRID_MAX) :: cdata
!!   integer :: isign
!!  USES
!! 
!!  AUTHOR
!!   C.M.Goringe/D.R.Bowler
!!  CREATION DATE
!!   12/09/95
!!  MODIFICATION HISTORY
!!   04/02/96 CMG
!!    Removed column_sorted
!!   01/06/2001 dave
!!    F90, portable datatypes, placed inside fft_module
!!    and removed include files
!!    Tidied up general use of modules
!!   06/06/2001 dave
!!    Added n_grid_? to use dimens statement
!!   06/06/2001 dave
!!    Changed call to rearrange_data
!!   24/06/2002 dave
!!    Changed calls to reflect gpfa use (default FFT package)
!!   08:29, 2003/09/05 dave
!!    Allocated trigs_x etc dynamically rather than statically
!!   2006/08/03 08:19 dave
!!    Passed in size of cdata, data
!!  SOURCE
!!
  subroutine fft3( data, cdata, size, isign ) 

    use dimens, ONLY: n_my_grid_points, n_grid_x, n_grid_y, n_grid_z
    use GenComms, ONLY: cq_abort, inode
    use maxima_module, ONLY: maxngrid
    use global_module, ONLY: area_SC
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl, type_int

    implicit none

    ! Passed variables
    integer :: size
    real(double), dimension(size) :: data
    complex(double_cplx), dimension(size):: cdata
    integer :: isign

    ! Local variables
    integer :: i, stat

    if(size>maxngrid) call cq_abort("Error in fft3 with sizes: ",size,maxngrid)
    if (isign==-1) then ! Forward transform (from domains to columns)
       do I = 1, n_my_grid_points
          cdata(I) = cmplx(data(I),0.0_double, double_cplx)
       end do
       call rearrange_data(cdata,maxngrid,packbx,psendbx,prbx,unpackbx)
       call fftx(cdata,maxngrid,x_columns_node(inode),n_grid_x,-1)
       call rearrange_data(cdata,maxngrid,packxy,psendxy,prxy,unpackxy)
       call ffty(cdata,maxngrid,y_columns_node(inode),n_grid_y,-1)
       call rearrange_data(cdata,maxngrid,packyz,psendyz,pryz,unpackyz)
       call fftz(cdata,maxngrid,z_columns_node(inode),n_grid_z,-1)
    else if (isign==1) then ! Reverse transform (from columns to domains)
       call fftz(cdata,maxngrid,z_columns_node(inode),n_grid_z,1)
       call rearrange_data(cdata,maxngrid,packzy,psendzy,przy,unpackzy)
       call ffty(cdata,maxngrid,y_columns_node(inode),n_grid_y,1)
       call rearrange_data(cdata,maxngrid,packyx,psendyx,pryx,unpackyx)
       call fftx(cdata,maxngrid,x_columns_node(inode),n_grid_x,1)
       call rearrange_data(cdata,maxngrid,packxb,psendxb,prxb,unpackxb)
       do I = 1, n_my_grid_points
          data(I) = real(cdata(I),double)
       end do
    else if (isign==0) then ! Initialise the fft tables
       ! Allocate space for trigonometric "twiddle" factors
       ! Note that the size of trigs_? is at least 2*(2^p+3^q+5^r), where n_grid_? factors 
       ! as 2^p*3^q*5^r.  If we set the size of trigs_? to 2*n_grid_? we'll be safe.
       allocate(trigs_x(2*n_grid_x),trigs_y(2*n_grid_y),trigs_z(2*n_grid_z),STAT=stat)
       if(stat/=0) call cq_abort("fft3: error allocating trigs ",stat)
       call reg_alloc_mem(area_SC,2*n_grid_x+2*n_grid_y+2*n_grid_z,type_dbl)
       call fftx(cdata,maxngrid,1,n_grid_x,0)
       call ffty(cdata,maxngrid,1,n_grid_y,0)
       call fftz(cdata,maxngrid,1,n_grid_z,0)
    end if
    return    
  end subroutine fft3
!!***

! -----------------------------------------------------------
! Subroutine rearrange_data
! -----------------------------------------------------------

!!****f* fft_module/rearrange_data *
!!
!!  NAME 
!!   rearrange_data
!!  USAGE
!! 
!!  PURPOSE
!!   Rearrange data in the manner described by the 
!!   pack,psnd,pr/rpr,unpack variables, created by set_fft_map.
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   C.M.Goringe/D.R.Bowler
!!  CREATION DATE
!!   06/03/95
!!  MODIFICATION HISTORY
!!   01/06/2001 dave
!!    F90, indentation, portable datatypes
!!   06/06/2001 dave
!!    Added use GenComms, and changed zmexch to exch
!!   2006/08/03 08:19 dave
!!    Passed in size of cdata
!!  SOURCE
!!
  subroutine rearrange_data( cdata, size, pack, psend, pr, unpack )

    use global_module, ONLY: numprocs, area_SC
    use maxima_module, ONLY: maxngrid
    use GenComms, ONLY: exchv, cq_abort
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl, type_int

    ! Passed variables
    integer :: size
    complex(double_cplx), dimension(size) :: cdata
    complex(double_cplx), allocatable, dimension(:) :: send
    complex(double_cplx), allocatable, dimension(:) :: recv

    integer :: pack(:), psend(:), pr(:), unpack(:)
    
    ! Local variables
    integer :: I, stat

    allocate(send(maxngrid), recv(maxngrid), STAT=stat)
    if(stat/=0) call cq_abort("Allocation error in rearrange_data: ",size,stat)
    call reg_alloc_mem(area_SC,4*maxngrid,type_dbl)
    ! pack data into sending format
    do I = 1, psend(numprocs+1)-1
       send(pack(I)) = cdata(I)
    end do

    ! call global exchange routine
    call exchv(send, psend, recv, pr, maxngrid)

    ! unpack data
    if(pr(numprocs+1)-1>size) call cq_abort("Error in rearrange_data: ",size,pr(numprocs+1)-1)
    do I = 1, pr(numprocs+1)-1
       cdata(I) = recv(unpack(I))
    end do
    deallocate(send, recv, STAT=stat)
    if(stat/=0) call cq_abort("Deallocation error in rearrange_data: ",size,stat)
    call reg_dealloc_mem(area_SC,4*maxngrid,type_dbl)
    return
  end subroutine rearrange_data
!!***

! -----------------------------------------------------------
! Subroutine fftx
! -----------------------------------------------------------

!!****f* fft_module/fftx *
!!
!!  NAME 
!!   fftx
!!  USAGE
!! 
!!  PURPOSE
!!   One dimensional fft of cols consecutive blocks of len data
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   C.M.Goringe
!!  CREATION DATE
!!   06/03/95
!!  MODIFICATION HISTORY
!!   01/06/2001 dave
!!    ROBODoc header, imported into fft_module
!!   24/06/2002 dave
!!    Now uses GPFA (which requires reals only...)
!!   08:31, 2003/09/05 dave
!!    Individual trigs used
!!   2006/08/03 08:19 dave
!!    Passed in size of cdata
!!  SOURCE
!!
  subroutine fftx( cdata, size, cols, len, isign )

    use numbers
    use global_module, ONLY: area_SC
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl, type_int
    use GenComms, ONLY: cq_abort

    implicit none

    ! Passed variables
    integer :: cols, len, isign, size
    complex(double_cplx), dimension(size) :: cdata
    real(double) :: scale

    ! Local variables
    integer :: I,j, stat
    real(double), allocatable, dimension(:) :: a,b

    allocate(a(len),b(len),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating space in fftx: ",len)
    call reg_alloc_mem(area_SC,2*len,type_dbl)
    scale = one/REAL(len,double)
    scale = SQRT(scale)
    if(isign==0) then
       call setgpfa(trigs_x,len)
    else
       do I = 0, cols-1
          do j=1,len
             a(j) = real(cdata(I*len+j),double)
             b(j) = aimag(cdata(I*len+j))
          end do
          call gpfa(a,b,trigs_x,1,1,len,1,isign)
          do j=1,len
             cdata(I*len+j) = scale*cmplx(a(j),b(j),double_cplx)
          end do
       end do
    endif
    deallocate(a,b,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating space in fftx: ",len)
    call reg_dealloc_mem(area_SC,2*len,type_dbl)
    return
  end subroutine fftx
!!***

! -----------------------------------------------------------
! Subroutine ffty
! -----------------------------------------------------------

!!****f* fft_module/ffty *
!!
!!  NAME 
!!   ffty
!!  USAGE
!! 
!!  PURPOSE
!!   One dimensional fft of cols consecutive blocks of len data
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   C.M.Goringe
!!  CREATION DATE
!!   06/03/95
!!  MODIFICATION HISTORY
!!   01/06/2001 dave
!!    ROBODoc header, imported into fft_module
!!   24/06/2002 dave
!!    Now uses GPFA (which requires reals only...)
!!   08:31, 2003/09/05 dave
!!    Individual trigs used
!!   2006/08/03 08:19 dave
!!    Passed in size of cdata
!!  SOURCE
!!
  subroutine ffty( cdata, size, cols, len, isign )

    use numbers
    use global_module, ONLY: area_SC
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl, type_int
    use GenComms, ONLY: cq_abort

    implicit none

    ! Passed variables
    integer :: cols, len, isign, size
    complex(double_cplx), dimension(size) :: cdata
    real(double) :: scale

    ! Local variables
    integer :: I,j, stat
    real(double), allocatable, dimension(:) :: a,b

    allocate(a(len),b(len),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating space in ffty: ",len)
    call reg_alloc_mem(area_SC,2*len,type_dbl)
    scale = one/REAL(len,double)
    scale = SQRT(scale)
    if(isign==0) then
       call setgpfa(trigs_y,len)
    else
       do I = 0, cols-1
          do j=1,len
             a(j) = real(cdata(I*len+j),double)
             b(j) = aimag(cdata(I*len+j))
          end do
          call gpfa(a,b,trigs_y,1,1,len,1,isign)
          do j=1,len
             cdata(I*len+j) = scale*cmplx(a(j),b(j),double_cplx)
          end do
       end do
    endif
    deallocate(a,b,STAT=stat)
    if(stat/=0) call cq_abort("Error allocating space in ffty: ",len)
    call reg_dealloc_mem(area_SC,2*len,type_dbl)
    return
  end subroutine ffty
!!***

! -----------------------------------------------------------
! Subroutine fftz
! -----------------------------------------------------------

!!****f* fft_module/fftz *
!!
!!  NAME 
!!   fftz
!!  USAGE
!! 
!!  PURPOSE
!!   One dimensional fft of cols consecutive blocks of len data
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   C.M.Goringe
!!  CREATION DATE
!!   06/03/95
!!  MODIFICATION HISTORY
!!   01/06/2001 dave
!!    ROBODoc header, imported into fft_module
!!   24/06/2002 dave
!!    Now uses GPFA (which requires reals only...)
!!   08:31, 2003/09/05 dave
!!    Individual trigs used
!!   2006/08/03 08:19 dave
!!    Passed in size of cdata
!!  SOURCE
!!
  subroutine fftz( cdata, size, cols, len, isign )

    use numbers
    use global_module, ONLY: area_SC
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl, type_int
    use GenComms, ONLY: cq_abort

    implicit none

    ! Passed variables
    integer :: cols, len, isign, size
    complex(double_cplx), dimension(size) :: cdata
    real(double) :: scale

    ! Local variables
    integer :: I,j, stat
    real(double), allocatable, dimension(:) :: a,b

    allocate(a(len),b(len),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating space in fftz: ",len)
    call reg_alloc_mem(area_SC,2*len,type_dbl)
    scale = one/REAL(len,double)
    scale = SQRT(scale)
    if(isign==0) then
       call setgpfa(trigs_z,len)
    else
       do I = 0, cols-1
          do j=1,len
             a(j) = real(cdata(I*len+j),double)
             b(j) = aimag(cdata(I*len+j))
          end do
          call gpfa(a,b,trigs_z,1,1,len,1,isign)
          do j=1,len
             cdata(I*len+j) = scale*cmplx(a(j),b(j),double_cplx)
          end do
       end do
    endif
    deallocate(a,b,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating space in fftz: ",len)
    call reg_dealloc_mem(area_SC,2*len,type_dbl)
    return
  end subroutine fftz
!!***

! -----------------------------------------------------------
! Subroutine set_fft_map
! -----------------------------------------------------------

!!****f* fft_module/set_fft_map *
!!
!!  NAME 
!!   set_fft_map
!!  USAGE
!! 
!!  PURPOSE
!!   The first sort of transformation is going from the domain structure
!!   to strips parallel to the x-axis : aka x-columns
!!   The grid as a whole is numbered (0..n_grid_x-1)
!!   this node has (i_domain_x-1..i_domain_x+n_grid_domain_x-2)
!!
!!   The x-columns will be distributed such that x,y is on node
!!   MOD( y + n_grid_y*z , numprocs ) + 1
!!
!!   and it will be, on that node, column number
!! ((y + z*n_grid_y ) / numprocs) + 1
!!
!!   We need to create three indexes. The first (packbx) takes all of our data,
!!   which occupies the first (n_grid_domain_x*n_grid_domain_y*n_grid_domain_z)
!!   elements of chden and points to where it should be stored in sendtemp,
!!   which is ordered such that the data to be sent to each recipient node 
!!   is contiguous. The second, psendbx, is the data needed by zmexch, 
!!   pointing at the start and finish in sendtemp of the data for each node, 
!!   and EITHER rprbx, the location on the remote node to put the data into 
!! (for SMAL style) OR prbx, the location on the local node to put data from 
!!   the remote node (for PVM style). The third is the unpacking equivalent to 
!!   the first, pointing to where the ith item in the new, x-column ordered, 
!!   chden is to be found in gettemp.
!!
!!   Inspection of the equivalencing in map.inc will show that most of the
!!   data for the inverse (xb) transformation is equivalent to the data for
!!   the bx transform
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   C.M.Goringe
!!  CREATION DATE
!!   02/04/96
!!  MODIFICATION HISTORY
!!   31/05/2001 dave
!!    Put into fft_module, added ROBODoc header, generalised data types
!!    and removed (I hope) the dependence on NODES = 2^n
!!   04/06/2001 dave
!!    Started to remove the IEOR lines
!!   06/06/2001 dave
!!    Bug fix - added n_grid_?, r_super_?_squared to dimens line
!!   07/06/2001 dave
!!    Changed alltoalls to exch and fixed copying problem
!!   07/06/2001 dave
!!    Changed stop to cq_abort and removed use mpi
!!   08/06/2001 dave
!!    Added my_barrier to use GenComms
!!   2005/11/14 08:12 dave
!!    Added lines to calculate and store reciprocal lattice vectors (for GGA)
!!    Also changed (slightly) calculation of hartree factor to use reciprocal lattice vectors
!!
!!    Note: this routine ASSUMES orthorhombic cells
!!   11:06, 14/11/2005 drb 
!!    Also changed check for zero (used to be r2.NE.0, now r2>very_small)
!!   15:54, 27/04/2007 drb 
!!    Changed recip_vector for consistency
!!  SOURCE
!!
  subroutine set_fft_map ( )

    use global_module, ONLY: numprocs, area_SC, iprint_SC
    use numbers
    use dimens, ONLY: n_my_grid_points, n_grid_x, n_grid_y, n_grid_z, &
         r_super_x_squared, r_super_y_squared, r_super_z_squared, &
         r_super_x, r_super_y, r_super_z
    use grid_index, ONLY: grid_point_x, grid_point_y, grid_point_z
    use GenComms, ONLY: gcopy_diff, exch, cq_abort, my_barrier, inode
    use maxima_module, ONLY: maxngrid
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl, type_int

    implicit none

    ! Local variables
    integer :: I, send_count(numprocs), count(numprocs), z, y, node_to, &
         lz, ly, x, lx, JNODE, xcol, get_count(numprocs), site, cnt, &
         col_num, N
    integer :: ierr
    integer, dimension(numprocs) :: sx_columns_node

    real(double) :: xmin2, ymin2, zmin2, r2
    real(double) :: xmin, ymin, zmin, cutoff, tmp_r2
    integer :: remote_n_my_grid_points
    integer, allocatable, dimension(:) :: remote_grid_point_x, remote_grid_point_y, remote_grid_point_z, &
         temp_kl

    ! Allocate variables
    allocate(prbx(numprocs+1), prxb(numprocs+1), prxy(numprocs+1), pryx(numprocs+1), &
         pryz(numprocs+1), przy(numprocs+1), STAT=ierr)
    if(ierr/=0) call cq_abort("Allocation error for pr in set_fft: ",numprocs+1,ierr)
    call reg_alloc_mem(area_SC,6*(numprocs+1),type_int)
    allocate(rprbx(numprocs), rprxb(numprocs), rprxy(numprocs), rpryx(numprocs), &
         rpryz(numprocs), rprzy(numprocs), STAT=ierr)
    if(ierr/=0) call cq_abort("Allocation error for rpr in set_fft: ",numprocs,ierr)
    call reg_alloc_mem(area_SC,6*numprocs,type_int)
    allocate(packbx(maxngrid),packxb(maxngrid),packxy(maxngrid),packyx(maxngrid), &
         packyz(maxngrid),packzy(maxngrid),STAT=ierr)
    if(ierr/=0) call cq_abort("Allocation error for pack in set_fft: ",maxngrid,ierr)
    call reg_alloc_mem(area_SC,6*maxngrid,type_int)
    allocate(x_columns_node(numprocs), y_columns_node(numprocs), z_columns_node(numprocs), STAT=ierr)
    if(ierr/=0) call cq_abort("Allocation error for columns in set_fft: ",numprocs,ierr)
    call reg_alloc_mem(area_SC,3*numprocs,type_int)
    allocate(hartree_factor(maxngrid), recip_vector(maxngrid,3), STAT=ierr)
    if(ierr/=0) call cq_abort("Allocation error for reciprocal factors in set_fft: ",maxngrid,ierr)
    call reg_alloc_mem(area_SC,4*maxngrid,type_dbl)
    allocate(remote_grid_point_x(maxngrid),remote_grid_point_y(maxngrid),remote_grid_point_z(maxngrid),&
         temp_kl(maxngrid),STAT=ierr)
    if(ierr/=0) call cq_abort("Allocation error for remote grids in set_fft: ",maxngrid,ierr)
    call reg_alloc_mem(area_SC,4*maxngrid,type_int)
    ! Map the equivalent pointers onto their arrays
    unpackxb => packbx
    unpackbx => packxb
    unpackxy => packyx
    unpackyx => packxy
    unpackyz => packzy
    unpackzy => packyz
    psendbx => prxb
    psendxb => prbx
    psendxy => pryx
    psendyx => prxy
    psendyz => przy
    psendzy => pryz

    ! first of all, count the number of data points to be SENT to each other
    ! node.
    call my_barrier
    do I = 1, numprocs
       send_count(I) = 0
       count(I) = 0
    end do
    do n = 1, n_my_grid_points
       z = grid_point_z(n)-1
       y = grid_point_y(n)-1
       node_to = MOD(y+n_grid_y*z,numprocs) + 1
       send_count(node_to) = send_count(node_to) + 1
    end do
    psendbx(1) = 1
    do I = 1, numprocs
       psendbx(I+1) = psendbx(I) + send_count(I)
       if(iprint_SC>3) write(io_lun,*) inode," proc sending bx to ",I,send_count(I)
    end do
    if(iprint_SC>3) write(io_lun,*) inode," bx total: ",psendbx(numprocs)+send_count(numprocs)


    ! Make the packbx
    do n = 1, n_my_grid_points
       z = grid_point_z(n)-1
       y = grid_point_y(n)-1
       node_to = MOD(y+n_grid_y*z,numprocs) + 1
       packbx(n) = count(node_to) + psendbx(node_to)
       count(node_to) = count(node_to) + 1
    end do

    ! We need to work out where on the remote node this
    ! data is going to need to be sent. 
    call exch(send_count, get_count, 1)

    ! prbx(I) is where we will unpack data from node I (for PVM)
    prbx(1) = 1
    do I = 1, numprocs
       prbx(I+1) = prbx(I) + get_count(I)  
    end do
    x_columns_node(inode) = (prbx(numprocs)+get_count(numprocs))/n_grid_x
    if(iprint_SC>3) write(io_lun,*) inode," xcols: ",prbx(numprocs)+get_count(numprocs),x_columns_node(inode)
    if (prbx(numprocs)+get_count(numprocs)>maxngrid+1) &
         call cq_abort("Grid overflow (bx) in FFT: ",prbx(numprocs)+get_count(numprocs), maxngrid)

    ! return this value to be rprbx (for SMAL)
    ! note that we can all return the inverse, prxb, because it is the same
    ! as psendbx
    call exch(prbx, rprbx, 1)
    call exch(prxb, rprxb, 1)
    ! This can probably be done better (?allgather?)
    do i=1,numprocs
       sx_columns_node(i) = x_columns_node(inode)
    enddo
    call exch(sx_columns_node, x_columns_node,1)
    if(iprint_SC>3) write(io_lun,*) inode," x columns: ",x_columns_node

    ! now we need to figure out the unpacking.
    ! for each node, JNODE, we loop over all points and identify the ones which
    ! are sending to us, and then identify 'site', where in our x-columns data
    ! that point will be stored.
    cnt = 0
    do jnode = 1, numprocs
       call gcopy_diff(remote_n_my_grid_points,n_my_grid_points,jnode)
       Call gcopy_diff(remote_grid_point_x, grid_point_x, jnode, maxngrid )
       Call gcopy_diff(remote_grid_point_y, grid_point_y, jnode, maxngrid )
       Call gcopy_diff(remote_grid_point_z, grid_point_z, jnode, maxngrid )
       do n = 1, remote_n_my_grid_points
          z = remote_grid_point_z(n)-1
          y = remote_grid_point_y(n)-1
          if (MOD(y+n_grid_y*z,numprocs).EQ.(inode-1)) then
             x = remote_grid_point_x(n)-1
             xcol = ((y+n_grid_y*z)/numprocs) + 1
             cnt = cnt + 1
             site = (xcol-1)*n_grid_x + x + 1
             unpackbx(site) = cnt
          end if
       end do
    end do

    ! now we need to do the same for the transform from x-columns to y-columns.
    ! since we are in x-column, item 'n' in chden is our column number
    ! (n/n_grid_x) + 1
    ! and is the 
    ! MOD(n,n_grid_x) + 1 th item in it.
    !
    ! The Ath column on this node is the (A-1)*numprocs+inode th column overall
    do I = 1, numprocs
       send_count(I) = 0
       count(I) = 0
    end do
    do n = 1, x_columns_node(inode)
       col_num = (n-1)*numprocs + inode - 1
       z = col_num / n_grid_y
       do x = 0, n_grid_x-1
          node_to = MOD(x+n_grid_x*z,numprocs) + 1
          send_count(node_to) = send_count(node_to) + 1
       end do
    end do
    psendxy(1) = 1
    do I = 1, numprocs
       psendxy(I+1) = psendxy(I) + send_count(I)
       if(iprint_SC>3) write(io_lun,*) inode," proc sending xy to ",I,send_count(I)
    end do
    if(iprint_SC>3) write(io_lun,*) inode," xy total: ",psendbx(numprocs)+send_count(numprocs)

    ! now we can make the packxy
    do n = 1, x_columns_node(inode)
       col_num = (n-1)*numprocs + inode - 1
       z = col_num / n_grid_y
       do x = 0, n_grid_x-1
          node_to = MOD(x+n_grid_x*z,numprocs) + 1
          packxy((n-1)*n_grid_x + x + 1) = count(node_to) + psendxy(node_to)
          count(node_to) = count(node_to) + 1
       end do
    end do

    ! now, we need to work out where on the remote node this
    ! data is going to need to be sent. Assume numprocs = 2**NDIMEN
    call exch(send_count, get_count, 1)

    ! prxy(I) is where we will unpack data from node I (for PVM)
    prxy(1) = 1
    do I = 1, numprocs
       prxy(I+1) = prxy(I) + get_count(I)  
    end do
    y_columns_node(inode) = (prxy(numprocs)+get_count(numprocs))/n_grid_y
    if (prxy(numprocs)+get_count(numprocs)>maxngrid+1) &
         call cq_abort("Grid overflow (xy) in FFT: ",prxy(numprocs)+get_count(numprocs), maxngrid)

    ! return this value to be rprxy (for SMAL)
    call exch(prxy, rprxy, 1)
    call exch(pryx, rpryx, 1)
    ! This can probably be done better
    do i=1,numprocs
       sx_columns_node(i) = y_columns_node(inode)
    end do
    call exch(sx_columns_node, y_columns_node,1)
    if(iprint_SC>3) write(io_lun,*) inode," y columns: ",y_columns_node

    ! now we need to figure out the unpacking.
    cnt = 0
    do jnode = 1, numprocs
       ! for each x-column stored on jnode
       do n = 1, x_columns_node(jnode)
          col_num = (n-1)*numprocs + jnode - 1
          y = MOD(col_num, n_grid_y)
          z = col_num / n_grid_y
          ! for each x point in that column
          do x = 0, n_grid_x-1
             node_to = MOD(x+n_grid_x*z,numprocs) + 1
             ! has it come to us?
             if (node_to.EQ.inode) then
                cnt = cnt + 1
                col_num = ((x+z*n_grid_x)/numprocs) + 1
                site = (col_num-1)*n_grid_y + y + 1
                unpackxy(site) = cnt
             end if
          end do
       end do
    end do

    ! finally we need to do the same for the transform from y-columns to 
    ! z-columns. In this case we also set up the scaling rules for the 
    ! hartree potential
    do I = 1, numprocs
       send_count(I) = 0
       count(I) = 0
    end do
    do n = 1, y_columns_node(inode)
       col_num = (n-1)*numprocs + inode - 1
       x = MOD(col_num,n_grid_x)
       do y = 0, n_grid_y-1
          node_to = MOD(y+n_grid_y*x,numprocs) + 1
          send_count(node_to) = send_count(node_to) + 1
       end do
    end do
    psendyz(1) = 1
    do I = 1, numprocs
       psendyz(I+1) = psendyz(I) + send_count(I)
       if(iprint_SC>3) write(io_lun,*) inode," proc sending yz to ",I,send_count(I)
    end do
    if(iprint_SC>3) write(io_lun,*) inode," yz total: ",psendbx(numprocs)+send_count(numprocs)

    ! now we can make the packyz
    do n = 1, y_columns_node(inode)
       col_num = (n-1)*numprocs + inode - 1
       x = MOD(col_num,n_grid_x)
       do y = 0, n_grid_y-1
          node_to = MOD(y+n_grid_y*x,numprocs) + 1
          packyz((n-1)*n_grid_y + y + 1) = count(node_to) + psendyz(node_to)
          count(node_to) = count(node_to) + 1
       end do
    end do

    ! now, we need to work out where on the remote node this
    ! data is going to need to be sent. Assume numprocs = 2**NDIMEN
    call exch(send_count, get_count, 1)

    ! pryz(I) is where we will unpack data from node I (for PVM)
    pryz(1) = 1
    do I = 1, numprocs
       pryz(I+1) = pryz(I) + get_count(I)  
    end do
    z_columns_node(inode) = (pryz(numprocs)+get_count(numprocs))/n_grid_z
    if(iprint_SC>3) write(io_lun,*) inode," z column: ",z_columns_node(inode)
    if (pryz(numprocs)+get_count(numprocs)>maxngrid+1) &
         call cq_abort("Grid overflow (yz) in FFT: ",pryz(numprocs)+get_count(numprocs), maxngrid)

    ! return this value to be rprxy (for SMAL)
    call exch(pryz, rpryz, 1)
    call exch(przy, rprzy, 1)
    ! This can probably be done better
    do i=1,numprocs
       sx_columns_node(i) = z_columns_node(inode)
    end do
    call exch(sx_columns_node, z_columns_node, 1)
    if(iprint_SC>3) write(io_lun,*) inode," z columns: ",z_columns_node

    ! Base kerker cutoff on shortest cell side and quarter
    if(r_super_x<r_super_y) then
       if(r_super_x<r_super_z) then
          cutoff = pi*n_grid_x/r_super_x
       else
          cutoff = pi*n_grid_z/r_super_z
       end if
    else
       if(r_super_y<r_super_z) then
          cutoff = pi*n_grid_y/r_super_y
       else
          cutoff = pi*n_grid_z/r_super_z
       end if
    end if
    ! This should be user-definable (currently assumes 0.5*longest vector)
    cutoff = 0.25_double*cutoff*cutoff
    size_kl = 0
    i0 = 0
    ! now we need to figure out the unpacking.
    cnt = 0
    do JNODE = 1, numprocs
       ! for each y-column stored on JNODE
       do n = 1, y_columns_node(JNODE)
          col_num = (n-1)*numprocs + JNODE - 1
          x = MOD(col_num,n_grid_x)
          !c xmin2 is distance to nearest edge of reciprocal box squared
          !xmin2 = real(min(x*x,(n_grid_x-x)*(n_grid_x-x)))
          !xmin2 = xmin2 / (r_super_x_squared)
          ! If we were looping over grid points in x, then I'd do a loop from 1 to ngridx/2+1 and set 
          ! x = (nx-1)/r_super_x and then a second loop from ngridx/2+2 to ngridx and set 
          ! x = ((nx-1)-ngridx)/r_super_x to get the Brillouin zone right; I'm assuming that the minus
          ! one is correctly done
          if((n_grid_x-x)<x) then
             xmin = real(x-n_grid_x,double)/r_super_x
          else
             xmin = real(x,double)/r_super_x
          end if
          xmin2 = xmin*xmin
          z = col_num / n_grid_x
          !zmin2 = real(min(z*z,(n_grid_z-z)*(n_grid_z-z)))
          !zmin2 = zmin2 / (r_super_z_squared)
          if((n_grid_z-z)<z) then
             zmin = real(z-n_grid_z,double)/r_super_z
          else
             zmin = real(z,double)/r_super_z
          end if
          zmin2 = zmin*zmin
          do y = 0, n_grid_y-1
             node_to = MOD(y+n_grid_y*x,numprocs) + 1
             if (node_to.EQ.inode) then
                cnt = cnt + 1
                col_num = ((y+x*n_grid_y)/numprocs) + 1
                site = (col_num-1)*n_grid_z + z + 1
                unpackyz(site) = cnt
                !ymin2 = real(min(y*y,(n_grid_y-y)*(n_grid_y-y)))
                !ymin2 =  ymin2 / (r_super_y_squared)
                if((n_grid_y-y)<y) then
                   ymin = real(y-n_grid_y,double)/r_super_y
                else
                   ymin = real(y,double)/r_super_y
                end if
                ymin2 = ymin*ymin
                recip_vector(site,1) = two*pi*xmin
                recip_vector(site,2) = two*pi*ymin
                recip_vector(site,3) = two*pi*zmin
                r2 = xmin2+ymin2+zmin2
                if (r2>very_small) then 
                   tmp_r2 = one/r2
                else
                   tmp_r2 = zero
                   i0 = site ! Store location of gamma point
                end if
                hartree_factor(site) = tmp_r2 
                if(tmp_r2<cutoff.AND.tmp_r2>very_small) then
                   size_kl = size_kl + 1
                   temp_kl(size_kl) = site
                end if
             end if
          end do
       end do
    end do
    allocate(kerker_list(size_kl),STAT=ierr)
    if(ierr/=0) call cq_abort("Error allocating kerker_list: ",size_kl,ierr)
    call reg_alloc_mem(area_SC,size_kl,type_int)
    do n=1,size_kl
       kerker_list(n) = temp_kl(n)
    end do
    deallocate(temp_kl,remote_grid_point_x,remote_grid_point_y,remote_grid_point_z,STAT=ierr)
    if(ierr/=0) call cq_abort("Deallocation error for remote grids in set_fft: ",maxngrid,ierr)
    call reg_dealloc_mem(area_SC,4*maxngrid,type_int)
    return
  end subroutine set_fft_map
!!***

end module fft_module
