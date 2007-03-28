! $Id: hartree.f90,v 1.6 2004/11/12 02:39:22 drb Exp $
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
!!  SOURCE
!!
module hartree_module

  implicit none

  ! RCS tag for object file identification
  character(len=80), save, private :: RCSid = "$Id: hartree.f90,v 1.6 2004/11/12 02:39:22 drb Exp $"

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
    !write(*,*) 'G=0 component is: ',chdenr(i0)
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

  subroutine kerker(resid,size,q0)

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

    !write(*,*) 'In kerker with q0: ',q0,size
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
    !write(*,*) 'Calling fft3'
    call fft3( resid, chdenr, size, -1 )
    !write(*,*) 'Called fft3'
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
    !write(*,*) 'Calling fft3'
    call fft3( resid, chdenr, size, +1 )
    !write(*,*) 'Called fft3'
    deallocate(chdenr)
    if(stat/=0) call cq_abort("Error deallocating chdenr in kerker: ",size,stat)
    call reg_dealloc_mem(area_SC,2*size,type_dbl)
    return
  end subroutine kerker

end module hartree_module
