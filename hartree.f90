! $Id$
! -----------------------------------------------------------
! Subroutine hartree
! -----------------------------------------------------------
! Code area 5: self-consistency
! -----------------------------------------------------------

!!***h* Conquest/hartree_module *
!!  NAME
!!   hartree_module
!!  CREATION DATE
!!   2004/10/05
!!  MODIFICATION HISTORY
!!   2004/10/29 drb
!!    Note: this was a subroutine, but was incorporated into a module to add kerker preconditioning
!!   2008/02/06 08:17 dave
!!    Changed for output to file not stdout
!!  SOURCE
!!
module hartree_module

  use global_module, ONLY: io_lun

  implicit none

  ! RCS tag for object file identification
  character(len=80), save, private :: RCSid = "$Id$"

!!***

contains


!!****f* hartree_module/hartree *
!!
!!  NAME 
!!   hartree
!!  USAGE
!! 
!!  PURPOSE
!!   Takes the charge density on the grid, normalised such that
!!   the sum over all points is equal to the number of electrons,
!!   and evaluates the hartree energy and potential. Assumes that
!!   set_fft_map has been called, and the fft tables have been initialised
!!  INPUTS
!!   real(double), dimension(N_GRID_MAX) :: chden - charge density
!!   real(double), dimension(N_GRID_MAX) :: potential - resulting potential
!!   real(double) :: energy - hartee energy due to charge
!!  USES
!! 
!!  AUTHOR
!!   C.M.Goringe
!!  CREATION DATE
!!   06/03/95
!!  MODIFICATION HISTORY
!!   31/05/2001 dave
!!    ROBODoc header, indented, added RCS Id tag
!!   01/06/2001 dave
!!    Fixed bug and added RCS Log tag as a test
!!    Added comments and tidied
!!   06/06/2001 dave
!!    Added dependence on fft_module
!!    Corrected dependence by adding z_columns node
!!   08/06/2001 dave
!!    Changed to use gsum from GenComms
!!   14:09, 29/07/2003 drb 
!!    Removed n_nodes from common use line
!!   11:38, 12/11/2004 dave 
!!    Removed grid_index use
!!  TODO
!!   Make absolutely sure about the following:
!!
!!   Originally, the subroutine started by scaling chden by 
!!   grid_point_volume (as the density is in e/au^3) at the
!!   start, kept factors of one_over_grid_point_volume through
!!   the routine (for dumr and dumi) and removed the factors from
!!   energy calculation, and then rescaled the charge at the end.
!!   I removed this, which adds a grid_point_volume multiple to
!!   the energy only. 31/05/2001 dave
!!  SOURCE
!!
  subroutine hartree( chden, potential, size, energy )

    use datatypes
    use dimens, ONLY: grid_point_volume, one_over_grid_point_volume, n_grid_z
    use numbers
    use fft_module, ONLY: fft3, hartree_factor, z_columns_node, i0
    use GenComms, ONLY: gsum,  inode, cq_abort
    use global_module, ONLY: area_SC
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    ! Passed variables
    integer :: size
    real(double), dimension(size), intent(in) :: chden
    real(double), dimension(size), intent(out) :: potential
    real(double) :: energy

    ! Local variables
    integer :: i, stat

    complex(double_cplx), allocatable, dimension(:) :: chdenr
    complex(double_cplx) :: t

    real(double) :: dumi, dumr, rp, ip
    ! refcoul is energy of two electrons seperated by one unit of distance.
    real(double), parameter :: refcoul = one
    ! harcon is the constant needed for energy and potential. It assumes that
    ! the hartree_factor correctly described the G vector at every point in
    ! reciprocal space.  MUST BE IN HARTREES
    real(double), parameter :: harcon = refcoul/pi

    allocate(chdenr(size), STAT=stat)
    if(stat/=0) call cq_abort("Error allocating chdenr in hartree: ",size,stat)
    call reg_alloc_mem(area_SC,size,type_dbl)
    ! N.B. The comment below is irrelevant, I think 01/06/2001 dave
    ! the density is in e/au^3, so convert to e/grid_point
    ! call scal( N_GRID_MAX, grid_point_volume, chden, 1 )
    call fft3( chden, chdenr, size, -1 )
    !write(io_lun,*) 'G=0 component is: ',chdenr(i0)
    energy = zero
    do i = 1, z_columns_node(inode)*n_grid_z
       !dumr = dble(chdenr(i))*hartree_factor(i)*harcon*  &
       !           one_over_grid_point_volume
       !dumi = dimag(chdenr(i)) * hartree_factor(i) * harcon *  &
       !           one_over_grid_point_volume 
       !energy = energy + dumr*dble(chdenr(i)) + dumi*dimag(chdenr(i))
       ! dumr = dble(chdenr(i))*hartree_factor(i)*harcon 
       ! dumi = dimag(chdenr(i))*hartree_factor(i)*harcon
       t = chdenr(i)*minus_i
       rp = real(chdenr(i),double)
       ! Alternative would be ip = aimag(chdenr(i))
       ip = real(t,double)
       dumr = rp*hartree_factor(i)*harcon 
       dumi = ip*hartree_factor(i)*harcon
       ! I could probably move the gpv lower 01/06/2001 dave
       energy = energy +  &
            grid_point_volume*dumr*rp +  &
            grid_point_volume*dumi*ip
       chdenr(i) = cmplx(dumr, dumi, double_cplx)
    end do
    call fft3( potential, chdenr, size, +1 )
    ! call scal( N_GRID_MAX, one_over_grid_point_volume, chden, 1 )
    call gsum(energy)
    energy = energy * half
    deallocate(chdenr, STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating chdenr in hartree: ",size,stat)
    call reg_dealloc_mem(area_SC,size,type_dbl)
    return
  end subroutine hartree
!!***

  subroutine kerker_obsolete(resid,size,q0)

    use datatypes
    use dimens, ONLY: grid_point_volume, one_over_grid_point_volume, n_grid_z
    use numbers
    use fft_module, ONLY: fft3, hartree_factor, z_columns_node, i0, kerker_list, size_kl
    use GenComms, ONLY: gsum, gmin, inode, cq_abort
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use global_module, ONLY: area_SC

    implicit none

    ! Passed variables
    integer :: size
    real(double), dimension(size) :: resid
    real(double) :: q0

    ! Local variables
    integer :: i, stat, ikl

    !complex(double_cplx), dimension(size) :: chdenr
    complex(double_cplx), allocatable, dimension(:) :: chdenr

    real(double) :: fac, facmin
    real(double) :: q02

    !write(io_lun,*) 'In kerker with q0: ',q0,size
    if(abs(q0)<1.0e-8_double) then
       return
    else
       q02 = q0*q0
    endif
    facmin = 1.0e8_double
    ! FFT residual
    allocate(chdenr(size), STAT=stat)
    if(stat/=0) call cq_abort("Error allocating chdenr in kerker: ",size,stat)
    call reg_alloc_mem(area_SC,2*size,type_dbl)
    !write(io_lun,*) 'Calling fft3'
    call fft3( resid, chdenr, size, -1 )
    !write(io_lun,*) 'Called fft3'
    do i = 1, z_columns_node(inode)*n_grid_z
       if(hartree_factor(i)>very_small) then
          fac = 1.0_double
          !if(hartree_factor(i)<cutoff) fac = (1.0_double + hartree_factor(i)*q02)
          fac = 1.0_double/(1.0_double + hartree_factor(i)*q02)
          facmin = min(fac,facmin)
       end if
       chdenr(i) = chdenr(i)*fac
    end do
    !!   kerker_list stores indices of points within reciprocal cutoff
    !do i=1,size_kl
    !   ikl = kerker_list(i)
    !   fac = 1.0_double/(1.0_double + hartree_factor(ikl)*q02)
    !   facmin = min(fac,facmin)
    !   chdenr(ikl) = chdenr(ikl)*fac
    !end do
    call gmin(facmin)
    ! Treat g=0 point separately
    if(i0>0) then
       chdenr(i0) = -chdenr(i0)*facmin
    end if
    !write(io_lun,*) 'Calling fft3'
    call fft3( resid, chdenr, size, +1 )
    !write(io_lun,*) 'Called fft3'
    deallocate(chdenr)
    if(stat/=0) call cq_abort("Error deallocating chdenr in kerker: ",size,stat)
    call reg_dealloc_mem(area_SC,2*size,type_dbl)
    return
  end subroutine kerker_obsolete




!!****f* hartree_module/kerker *
!!
!!  NAME 
!!   kerker
!!  USAGE
!! 
!!  PURPOSE 
!!   take charge density residue from real space grid, fourier
!!   transform and then add kerker preconditioning factor
!!   correspodingly and fourier transform back
!!  INPUTS
!!   real(double), dimension(size) :: resid      - residue and results in Kerker preconditioned residue
!!   integer                       :: size       - maximum size of resid and resid_cov array
!!   real(double)                  :: q0         - Kerker preconditioning factor
!!  USES
!! 
!!  AUTHOR
!!   David Bowler, Lianheng Tong
!!  CREATION DATE
!!   2010/07/30
!!  MODIFICATION HISTORY
!!   2010/07/30  Lianheng Tong  
!!     Based on kerker_obsolete written by dave. 
!!  TODO
!!
!!  SOURCE
!!  
  subroutine kerker (resid, size, q0)
    
    use datatypes
    use dimens, ONLY: grid_point_volume, one_over_grid_point_volume, n_grid_z
    use numbers, ONLY: very_small, zero, one
    use fft_module, ONLY: fft3, hartree_factor, z_columns_node, i0, kerker_list, size_kl
    use GenComms, ONLY: gsum, inode, cq_abort
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use global_module, ONLY: area_SC

    implicit none
    
    ! Passed variables
    integer :: size
    real(double), dimension(size) :: resid
    real(double), intent(in) :: q0

    ! Local variables
    integer :: i, stat
    complex(double_cplx), allocatable, dimension(:) :: FR_kerker
    real(double) :: fac, q02
    
    if (abs(q0) < very_small) then
       return
    else
       q02 = q0*q0
    endif
    ! FFT residue
    allocate (FR_kerker(size), STAT=stat)
    if (stat /= 0) call cq_abort ("hartree_module/kerker: Error allocating FR_kerker: ", size, stat)
    call reg_alloc_mem (area_SC, size, type_dbl)
    call fft3 (resid, FR_kerker, size, -1)
    do i = 1, z_columns_node(inode)*n_grid_z
       ! hartree_factor(q) = 1/q**2 for q /= 0, and hartree_factor(q) = 0 for q = 0, calculated in fft module
       ! excluding q=0 point, treating it separately
       if (hartree_factor(i) > very_small) then 
          fac = one / (one + hartree_factor(i)*q02)
          FR_kerker(i) = FR_kerker(i)*fac
       end if
    end do
    ! do q=0 point separately (if q=0 is on one of the grid point)
    ! note that if q=0 point is not on the discrete reciporical grid for FFT, then i0=0
    ! q=0 is only on one of processor node need to make sure we are doing the correct point
    if ((i0 > 0) .AND. (hartree_factor(i0) <= very_small)) then
       FR_kerker(i0) = zero
    end if
    ! FFT back
    call fft3 (resid, FR_kerker, size, +1)
    ! deallocate array
    deallocate (FR_kerker, STAT=stat)
    if (stat /= 0) call cq_abort ("hartree_module/kerker: Error deallocating FR_kerker: ", size, stat)
    call reg_dealloc_mem (area_SC, size, type_dbl)

    return

  end subroutine kerker


!!****f* hartree_module/wdmetric *
!!
!!  NAME 
!!   wdmetric
!!  USAGE
!! 
!!  PURPOSE 
!!   take charge density residue from real space grid, fourier
!!   transform and then add wave-dependent-metric factor correspodingly and
!!   fourier transform back to an different array
!!  INPUTS
!!   real(double), dimension(size) :: resid      - residue input
!!   real(double), dimension(size) :: resid_cov  - results covariant residue for wave-dependent-metric
!!   integer                       :: size       - maximum size of resid and resid_cov array
!!   real(double)                  :: q1         - wave-dependent-metric factor
!!  USES
!! 
!!  AUTHOR
!!   Lianheng Tong
!!  CREATION DATE
!!   2010/07/30
!!  MODIFICATION HISTORY
!!
!!  TODO
!!
!!  SOURCE
!!  
  subroutine wdmetric (resid, resid_cov, size, q1)
    
    use datatypes
    use dimens, ONLY: grid_point_volume, one_over_grid_point_volume, n_grid_z
    use numbers, ONLY: very_small, zero, one
    use fft_module, ONLY: fft3, hartree_factor, z_columns_node, i0, kerker_list, size_kl
    use GenComms, ONLY: gsum, gmax, inode, cq_abort
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use global_module, ONLY: area_SC

    implicit none
    
    ! Passed variables
    integer :: size
    real(double), dimension(size), intent(in) :: resid
    real(double), dimension(size), intent(out) :: resid_cov
    real(double), intent(in) :: q1

    ! Local variables
    integer :: i, stat
    complex(double_cplx), allocatable, dimension(:) :: FR_wdmetric
    real(double) :: fac, facmax, q12
    
    if (abs (q1) < very_small) then
       ! copy resid to resid_cov 
       resid_cov = resid
       return
    else
       q12 = q1*q1
    end if
    facmax = zero
    ! FFT residue
    allocate (FR_wdmetric(size), STAT=stat)
    if (stat /= 0) call cq_abort ("hartree_module/wdmetric: Error allocating FR_wdmetric: ", size, stat)
    call reg_alloc_mem (area_SC, size, type_dbl)
    call fft3 (resid, FR_wdmetric, size, -1)
    do i = 1, z_columns_node(inode)*n_grid_z
       ! hartree_factor(q) = 1/q**2 for q /= 0, and hartree_factor(q) = 0 for q = 0, calculated in fft module
       ! excluding q=0 point, treating it separately
       if (hartree_factor(i) > very_small) then 
          fac = one + hartree_factor(i)*q12
          facmax = max (fac, facmax)
          FR_wdmetric(i) = FR_wdmetric(i)*fac
       end if
    end do
    ! find the global maximum for fac
    call gmax (facmax)
    ! do q=0 point separately (if q=0 is on one of the grid point)
    ! note that if q=0 point is not on the discrete reciporical grid for FFT, then i0=0
    ! q=0 is only on one of processor node need to make sure we are doing the correct point
    if ((i0 > 0) .AND. (hartree_factor(i0) <= very_small)) then
       FR_wdmetric(i0) = FR_wdmetric(i0)*facmax
    end if
    ! FFT back
    call fft3 (resid_cov, FR_wdmetric, size, +1)
    ! deallocate arrays
    deallocate (FR_wdmetric, STAT=stat)
    if (stat /= 0) call cq_abort ("hartree_module/wdmetric: Error deallocating FR_wdmetric: ", size, stat)
    call reg_dealloc_mem (area_SC, size, type_dbl)

    return

  end subroutine wdmetric



!!****f* hartree_module/kerker_and_wdmetric *
!!
!!  NAME 
!!   kerker_and_wdmetric
!!  USAGE
!! 
!!  PURPOSE 
!!   take charge density residue from real space grid, fourier
!!   transform and then add kerker preconditioning and wave-dependent-metric
!!   factor correspodingly (stored to a different array) and fourier
!!   transform back
!!  INPUTS
!!   real(double), dimension(size) :: resid      - residue and results in Kerker preconditioned residue
!!   real(double), dimension(size) :: resid_cov  - results covariant residue for wave-dependent-metric
!!   integer                       :: size       - maximum size of resid and resid_cov array
!!   real(double)                  :: q0         - Kerker preconditioning factor
!!   real(double)                  :: q1         - wave-dependent-metric factor
!!  USES
!! 
!!  AUTHOR
!!   Lianheng Tong
!!  CREATION DATE
!!   2010/07/30
!!  MODIFICATION HISTORY
!!
!!  TODO
!!
!!  SOURCE
!!  
  subroutine kerker_and_wdmetric (resid, resid_cov, size, q0, q1)
    
    use datatypes
    use dimens, ONLY: grid_point_volume, one_over_grid_point_volume, n_grid_z
    use numbers, ONLY: very_small, zero, one
    use fft_module, ONLY: fft3, hartree_factor, z_columns_node, i0, kerker_list, size_kl
    use GenComms, ONLY: gsum, gmax, inode, cq_abort
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use global_module, ONLY: area_SC

    implicit none
    
    ! Passed variables
    integer :: size
    real(double), dimension(size) :: resid
    real(double), dimension(size), intent(out) :: resid_cov
    real(double), intent(in) :: q0
    real(double), intent(in) :: q1

    ! Local variables
    integer :: i, stat
    complex(double_cplx), allocatable, dimension(:) :: FR_kerker, FR_wdmetric
    real(double) :: fac, fac2, facmax, q02, q12
    
    if (abs(q0) < very_small) then
       return
    else
       q02 = q0*q0
    endif
    if (abs (q1) < very_small) then
       ! copy resid to resid_cov 
       resid_cov = resid
       return
    else
       q12 = q1*q1
    end if
    facmax = zero
    ! FFT residue
    allocate (FR_kerker(size), FR_wdmetric(size), STAT=stat)
    if (stat /= 0) call cq_abort ("hartree_module/kerker_and_wdmetric: Error allocating FR_kerker, FR_wdmetric: ", size, stat)
    call reg_alloc_mem (area_SC, 2*size, type_dbl)
    call fft3 (resid, FR_kerker, size, -1)
    FR_wdmetric = FR_kerker
    do i = 1, z_columns_node(inode)*n_grid_z
       ! hartree_factor(q) = 1/q**2 for q /= 0, and hartree_factor(q) = 0 for q = 0, calculated in fft module
       ! excluding q=0 point, treating it separately
       if (hartree_factor(i) > very_small) then 
          ! Kerker factor
          fac = one / (one + hartree_factor(i)*q02)
          FR_kerker(i) = FR_kerker(i)*fac
          ! wave-dependent-metric factor
          fac2 = one + hartree_factor(i)*q12
          facmax = max (fac2, facmax)
          FR_wdmetric(i) = FR_wdmetric(i)*fac2
       end if
    end do
    ! find the global maximum for fac
    call gmax (facmax)
    ! do q=0 point separately (if q=0 is on one of the grid point)
    ! note that if q=0 point is not on the discrete reciporical grid for FFT, then i0=0
    ! q=0 is only on one of processor node need to make sure we are doing the correct point
    if ((i0 > 0) .AND. (hartree_factor(i0) <= very_small)) then
       FR_kerker(i0) = zero
       FR_wdmetric(i0) = FR_wdmetric(i0)*facmax
    end if
    ! FFT back
    call fft3 (resid, FR_kerker, size, +1)
    call fft3 (resid_cov, FR_wdmetric, size, +1)
    ! deallocate arrays
    deallocate (FR_kerker, FR_wdmetric, STAT=stat)
    if (stat /= 0) call cq_abort ("hartree_module/kerker_and_wdmetric: Error deallocating FR_kerker, FR_wdmetric: ", &
         size, stat)
    call reg_dealloc_mem (area_SC, 2*size, type_dbl)

    return

  end subroutine kerker_and_wdmetric

  



end module hartree_module
