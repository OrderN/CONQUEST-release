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
!!    Note: this was a subroutine, but was incorporated into a module
!!    to add kerker preconditioning
!!   2008/02/06 08:17 dave
!!    Changed for output to file not stdout
!!   2015/05/01 09:15 dave and sym
!!    Implementing stress (added Hartree stress as module variable)
!!   2019/04/08 zamaan
!!    Added off diagonal elements of hartree stress 
!!  SOURCE
!!
module hartree_module

  use global_module, only: io_lun

  use datatypes
  
  implicit none

  ! We have this here since hartree_module is used by force_module, so we can't have it here
  real(double), dimension(3,3) :: Hartree_stress

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
  !!   2016/01/21 08:26 dave
  !!    Notes should have been added during first implementation of stress
  !!    We now calculate stress from change in G with cell size (hartree stress)
  !!    as well as (optionally) stress for local pseudopotential coming from the
  !!    FFT of charge to make the potential (this is second stress)
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
  subroutine hartree(chden, potential, size, energy, chden2, stress2)

    use datatypes
    use numbers    
    use dimens,        only: grid_point_volume, &
                             one_over_grid_point_volume, n_grid_z
    use fft_module,    only: fft3, hartree_factor, z_columns_node, i0, recip_vector
    use GenComms,      only: gsum,  inode, cq_abort
    use global_module, only: area_SC, flag_full_stress, flag_stress
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    ! Passed variables
    integer      :: size
    real(double) :: energy
    real(double), dimension(size), intent(in)  :: chden
    real(double), dimension(size), intent(out) :: potential
    real(double), dimension(3,3) :: chden_str_r, chden_str_i
    real(double), dimension(size), OPTIONAL :: chden2
    real(double), dimension(3,3), OPTIONAL :: stress2

    ! Local variables
    integer :: i, stat, dir1, dir2

    complex(double_cplx), allocatable, dimension(:) :: chdenr
    complex(double_cplx), allocatable, dimension(:) :: str_chdenr

    real(double) :: dumi, dumr, rp, ip, rp2, ip2, rv2
    ! refcoul is energy of two electrons seperated by one unit of distance.
    real(double), parameter :: refcoul = one
    ! harcon is the constant needed for energy and potential. It assumes that
    ! the hartree_factor correctly described the G vector at every point in
    ! reciprocal space.  MUST BE IN HARTREES
    real(double), parameter :: harcon = refcoul/pi
    logical :: second_stress

    second_stress = .false.
    if(PRESENT(chden2).AND.PRESENT(stress2)) then
       second_stress = .true.
       allocate(str_chdenr(size), STAT=stat)
       call fft3(chden2, str_chdenr, size, -1)
       stress2 = zero
    end if    
    allocate(chdenr(size), STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error allocating chdenr in hartree: ", size, stat)
    call reg_alloc_mem(area_SC, size, type_dbl)
    ! N.B. The comment below is irrelevant, I think 01/06/2001 dave
    ! the density is in e/au^3, so convert to e/grid_point
    ! call scal( N_GRID_MAX, grid_point_volume, chden, 1 )
    call fft3(chden, chdenr, size, -1)
    !write(io_lun,*) 'G=0 component is: ',chdenr(i0)
    energy = zero
    Hartree_stress = zero
    if(second_stress) stress2 = zero
    do i = 1, z_columns_node(inode)*n_grid_z
       rp = real(chdenr(i), double)
       ! Originally, we had t=chdenr(i) * minus_i and ip = real(t,double) 
       ip = aimag(chdenr(i))
       dumr = rp * hartree_factor(i) 
       dumi = ip * hartree_factor(i) 
       energy = energy + dumr * rp + dumi * ip
       chdenr(i) = cmplx(dumr, dumi, double_cplx)
       if (flag_stress) then
         if(second_stress) then
            rp2 = real(str_chdenr(i),double)
            ip2 = aimag(str_chdenr(i))
         end if
       end if
       ! SYM 2014/09/08 15:41 Hartree stress calculate and accumulate
       if (flag_stress) then
         do dir1 = 1,3
            if (flag_full_stress) then
               do dir2 = 1,3
                 rv2 = recip_vector(i,dir1)*recip_vector(i,dir2)
                 chden_str_r(dir1,dir2) = dumr * rv2 * hartree_factor(i)
                 chden_str_i(dir1,dir2) = dumi * rv2 * hartree_factor(i)
                 Hartree_stress(dir1,dir2) = Hartree_stress(dir1,dir2) + &
                   rp*chden_str_r(dir1,dir2) + ip*chden_str_i(dir1,dir2)
                 if (second_stress) stress2(dir1,dir2) = stress2(dir1,dir2) + &
                   rp2*chden_str_r(dir1,dir2) + ip2*chden_str_i(dir1,dir2)
               end do
            else
               rv2 = recip_vector(i,dir1)*recip_vector(i,dir1)
               chden_str_r(dir1,dir1) = dumr * rv2 * hartree_factor(i)
               chden_str_i(dir1,dir1) = dumi * rv2 * hartree_factor(i)
               Hartree_stress(dir1,dir1) = Hartree_stress(dir1,dir1) + &
                 rp*chden_str_r(dir1,dir1) + ip*chden_str_i(dir1,dir1)
               if (second_stress) stress2(dir1,dir1) = stress2(dir1,dir1) + &
                 rp2*chden_str_r(dir1,dir1) + ip2*chden_str_i(dir1,dir1)
            end if
         end do
       end if
    end do
    call fft3(potential, chdenr, size, +1)
    ! Sum over processes
    call gsum(energy)
    if (flag_stress) call gsum(Hartree_stress,3,3)
    ! Scale
    potential = potential*harcon
    energy = energy * grid_point_volume*half* harcon 
    if (flag_stress) then
      Hartree_stress(1:3,1:3) = Hartree_stress(1:3,1:3) * harcon * &
                                grid_point_volume/(two*pi*two*pi)
      if(second_stress) then
         call gsum(stress2,3,3)
         stress2(1:3,1:3) = stress2(1:3,1:3) * two * harcon * &
                            grid_point_volume/(two*pi*two*pi)
      end if
    end if
    deallocate(chdenr, STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error deallocating chdenr in hartree: ", size, stat)
    call reg_dealloc_mem(area_SC, size, type_dbl)
    return
  end subroutine hartree
  !!***

  !!****f* hartree_module/get_hartree_energy *
  !!
  !!  NAME 
  !!   get_hartree_energy
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Takes the charge density on the grid, normalised such that
  !!   the sum over all points is equal to the number of electrons,
  !!   and evaluates the hartree energy
  !!  INPUTS
  !!   real(double), dimension(N_GRID_MAX) :: chden - charge density
  !!   real(double) :: energy - hartee energy due to charge
  !!  USES
  !! 
  !!  AUTHOR
  !!   D. R. Bowler
  !!  CREATION DATE
  !!   2021/07/22
  !!  MODIFICATION HISTORY
  !!
  !!  SOURCE
  !!
  subroutine get_hartree_energy(chden, size, energy)

    use datatypes
    use numbers    
    use dimens,        only: grid_point_volume, &
                             one_over_grid_point_volume, n_grid_z
    use fft_module,    only: fft3, hartree_factor, z_columns_node, i0, recip_vector
    use GenComms,      only: gsum,  inode, cq_abort
    use global_module, only: area_SC, flag_full_stress, flag_stress
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    ! Passed variables
    integer      :: size
    real(double) :: energy
    real(double), dimension(size), intent(in)  :: chden

    ! Local variables
    integer :: i, stat, dir1, dir2

    complex(double_cplx), allocatable, dimension(:) :: chdenr

    real(double) :: dumi, dumr, rp, ip, rp2, ip2, rv2
    ! refcoul is energy of two electrons seperated by one unit of distance.
    real(double), parameter :: refcoul = one
    ! harcon is the constant needed for energy and potential. It assumes that
    ! the hartree_factor correctly described the G vector at every point in
    ! reciprocal space.  MUST BE IN HARTREES
    real(double), parameter :: harcon = refcoul/pi
    logical :: second_stress

    allocate(chdenr(size), STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error allocating chdenr in get_hartree_energy: ", size, stat)
    call reg_alloc_mem(area_SC, 2*size, type_dbl)
    ! FFT
    call fft3(chden, chdenr, size, -1)
    energy = zero
    do i = 1, z_columns_node(inode)*n_grid_z
       ! Energy is sum over n(G)^2/G^2
       rp = real(chdenr(i), double)
       ip = aimag(chdenr(i))
       energy = energy + (rp * rp + ip * ip) * hartree_factor(i)
    end do
    ! Sum over processes
    call gsum(energy)
    energy = energy * grid_point_volume * half * harcon 
    deallocate(chdenr, STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error deallocating chdenr in get_hartree_energy: ", size, stat)
    call reg_dealloc_mem(area_SC, 2*size, type_dbl)
    return
  end subroutine get_hartree_energy
  !!***

  subroutine kerker_obsolete(resid,size,q0)

    use datatypes
    use dimens,        only: grid_point_volume,                        &
                             one_over_grid_point_volume, n_grid_z
    use numbers
    use fft_module,    only: fft3, hartree_factor, z_columns_node, i0, &
                             kerker_list, size_kl
    use GenComms,      only: gsum, gmin, inode, cq_abort
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use global_module, only: area_SC

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
    if(abs(q0)<RD_ERR) then
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
    call fft3(resid, chdenr, size, -1)
    !write(io_lun,*) 'Called fft3'
    do i = 1, z_columns_node(inode)*n_grid_z
       if(hartree_factor(i)>RD_ERR) then
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
    call fft3(resid, chdenr, size, +1)
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
  !!   real(double), dimension(size) :: resid  - residue and results in Kerker
  !!                                             preconditioned  residue
  !!   integer                       :: size   - maximum size of
  !!                                             resid and resid_cov array
  !!   real(double)                  :: q0     - Kerker preconditioning factor
  !!  USES
  !! 
  !!  AUTHOR
  !!   David Bowler, Lianheng Tong
  !!  CREATION DATE
  !!   2010/07/30
  !!  MODIFICATION HISTORY
  !!   2010/07/30  Lianheng Tong  
  !!     Based on kerker_obsolete written by dave. 
  !!   2011/07/28 14:20 dave and umberto
  !!     Fix problem with AND statement referencing hartree_factor
  !!     with index 0
  !!   2023/06/02 16:04 dave
  !!     Remove scaling for q=0 point (introduces errors)
  !!  TODO
  !!
  !!  SOURCE
  !!  
  subroutine kerker(resid, size, q0)
    
    use datatypes
    use dimens,        only: grid_point_volume,                        &
                             one_over_grid_point_volume, n_grid_z
    use numbers,       only: RD_ERR, zero, one
    use fft_module,    only: fft3, hartree_factor, z_columns_node, i0, &
                             kerker_list, size_kl
    use GenComms,      only: gsum, inode, cq_abort
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use global_module, only: area_SC

    implicit none
    
    ! Passed variables
    integer,      intent(in) :: size
    real(double), intent(in) :: q0
    real(double), dimension(size) :: resid

    ! Local variables
    integer      :: i, stat
    real(double) :: fac, q02
    complex(double_cplx), allocatable, dimension(:) :: FR_kerker
    
    if (abs(q0) < RD_ERR) then
       return
    else
       q02 = q0*q0
    endif
    ! FFT residue
    allocate(FR_kerker(size), STAT=stat)
    if (stat /= 0) &
         call cq_abort("hartree_module/kerker: Error allocating FR_kerker: ", &
                       size, stat)
    call reg_alloc_mem(area_SC, size, type_dbl)
    call fft3(resid, FR_kerker, size, -1)
    do i = 1, z_columns_node(inode)*n_grid_z
       ! hartree_factor(q) = 1/q**2 for q /= 0, and hartree_factor(q)
       ! = 0 for q = 0, calculated in fft module
       if (hartree_factor(i) > RD_ERR) then 
          fac = one / (one + hartree_factor(i)*q02)
          FR_kerker(i) = FR_kerker(i)*fac
       end if
    end do
    ! do nothing for q=0 point (introduces an error!)
    ! FFT back
    call fft3(resid, FR_kerker, size, +1)
    ! deallocate array
    deallocate(FR_kerker, STAT=stat)
    if (stat /= 0) &
         call cq_abort("hartree_module/kerker: Error deallocating FR_kerker: ", &
                       size, stat)
    call reg_dealloc_mem(area_SC, size, type_dbl)

    return

  end subroutine kerker
  !*****


  !!****f* hartree_module/wdmetric *
  !!
  !!  NAME 
  !!   wdmetric
  !!  USAGE
  !! 
  !!  PURPOSE 
  !!   take charge density residue from real space grid, fourier
  !!   transform and then add wave-dependent-metric factor
  !!   correspodingly and fourier transform back to an different array
  !!  INPUTS
  !!   real(double), dimension(size) :: resid      - residue input
  !!   real(double), dimension(size) :: resid_cov  - results covariant
  !!                                                 residue for wave-
  !!                                                 dependent-metric
  !!   integer                       :: size       - maximum size of
  !!                                                 resid and resid_cov array
  !!   real(double)                  :: q1         - wave-dependent-metric factor
  !!  USES
  !! 
  !!  AUTHOR
  !!   Lianheng Tong
  !!  CREATION DATE
  !!   2010/07/30
  !!  MODIFICATION HISTORY
  !!   2011/07/28 14:20 dave and umberto
  !!     Fix problem with AND statement referencing hartree_factor
  !!     with index 0
  !!   2023/06/02 16:04 dave
  !!     Remove scaling for q=0 point (introduces errors)
  !!  TODO
  !!
  !!  SOURCE
  !!  
  subroutine wdmetric(resid, resid_cov, size, q1)
    
    use datatypes
    use dimens,        only: grid_point_volume,                        &
                             one_over_grid_point_volume, n_grid_z
    use numbers,       only: RD_ERR, zero, one
    use fft_module,    only: fft3, hartree_factor, z_columns_node, i0, &
                             kerker_list, size_kl
    use GenComms,      only: gsum, gmax, inode, cq_abort
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use global_module, only: area_SC

    implicit none
    
    ! Passed variables
    integer,      intent(in) :: size
    real(double), intent(in) :: q1
    real(double), dimension(size), intent(in)  :: resid
    real(double), dimension(size), intent(out) :: resid_cov

    ! Local variables
    integer      :: i, stat
    real(double) :: fac, q12
    complex(double_cplx), allocatable, dimension(:) :: FR_wdmetric
    
    if (abs (q1) < RD_ERR) then
       ! copy resid to resid_cov 
       resid_cov = resid
       return
    else
       q12 = q1*q1
    end if
    ! FFT residue
    allocate(FR_wdmetric(size), STAT=stat)
    if (stat /= 0) &
         call cq_abort("hartree_module/wdmetric: Error allocating FR_wdmetric: ",&
                       size, stat)
    call reg_alloc_mem(area_SC, size, type_dbl)
    call fft3(resid, FR_wdmetric, size, -1)
    do i = 1, z_columns_node(inode)*n_grid_z
       ! hartree_factor(q) = 1/q**2 for q /= 0, and hartree_factor(q)
       ! = 0 for q = 0, calculated in fft module excluding q=0 point,
       if (hartree_factor(i) > RD_ERR) then 
          fac = one + hartree_factor(i)*q12
          FR_wdmetric(i) = FR_wdmetric(i)*fac
       end if
    end do
    ! do nothing for q=0 point (nothing needed)
    ! FFT back
    call fft3(resid_cov, FR_wdmetric, size, +1)
    ! deallocate arrays
    deallocate(FR_wdmetric, STAT=stat)
    if (stat /= 0) &
         call cq_abort("hartree_module/wdmetric: Error deallocating FR_wdmetric: ",&
                       size, stat)
    call reg_dealloc_mem(area_SC, size, type_dbl)

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
  !!   real(double), dimension(size) :: resid      - residue and results
  !!                                                 in Kerker preconditioned
  !!                                                 residue
  !!   real(double), dimension(size) :: resid_cov  - results covariant residue
  !!                                                 for wave-dependent-metric
  !!   integer                       :: size       - maximum size of resid and
  !!                                                 resid_cov array
  !!   real(double)                  :: q0         - Kerker preconditioning factor
  !!   real(double)                  :: q1         - wave-dependent-metric factor
  !!  USES
  !! 
  !!  AUTHOR
  !!   Lianheng Tong
  !!  CREATION DATE
  !!   2010/07/30
  !!  MODIFICATION HISTORY
  !!   2011/07/28 14:20 dave and umberto
  !!     Fix problem with AND statement referencing hartree_factor
  !!     with index 0
  !!   2023/06/02 16:04 dave
  !!     Remove scaling for q=0 point (introduces errors)
  !!  TODO
  !!
  !!  SOURCE
  !!  
  subroutine kerker_and_wdmetric(resid, resid_cov, size, q0, q1)
    
    use datatypes
    use dimens,        only: grid_point_volume,                        &
                             one_over_grid_point_volume, n_grid_z
    use numbers,       only: RD_ERR, zero, one
    use fft_module,    only: fft3, hartree_factor, z_columns_node, i0, &
                             kerker_list, size_kl
    use GenComms,      only: gsum, gmax, inode, cq_abort
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use global_module, only: area_SC

    implicit none
    
    ! Passed variables
    integer,      intent(in) :: size
    real(double), intent(in) :: q0
    real(double), intent(in) :: q1
    real(double), dimension(size) :: resid
    real(double), dimension(size), intent(out) :: resid_cov

    ! Local variables
    integer      :: i, stat
    real(double) :: fac, fac2, q02, q12
    complex(double_cplx), allocatable, dimension(:) :: FR_kerker, FR_wdmetric
    
    q02 = q0*q0
    q12 = q1*q1
    ! FFT residue
    allocate(FR_kerker(size), FR_wdmetric(size), STAT=stat)
    if (stat /= 0) &
         call cq_abort("hartree_module/kerker_and_wdmetric: Error &
                        &allocating FR_kerker, FR_wdmetric: ",    &
                       size, stat)
    call reg_alloc_mem(area_SC, 2*size, type_dbl)
    call fft3(resid, FR_kerker, size, -1)
    FR_wdmetric = FR_kerker
    do i = 1, z_columns_node(inode)*n_grid_z
       ! hartree_factor(q) = 1/q**2 for q /= 0, and hartree_factor(q)
       ! = 0 for q = 0, calculated in fft module excluding q=0 point,
       if (hartree_factor(i) > RD_ERR) then 
          ! Kerker factor
          fac = one / (one + hartree_factor(i)*q02)
          FR_kerker(i) = FR_kerker(i)*fac
          ! wave-dependent-metric factor
          fac2 = one + hartree_factor(i)*q12
          FR_wdmetric(i) = FR_wdmetric(i)*fac2
       end if
    end do
    ! do nothing for q=0 point (nothing needed)
    ! FFT back
    call fft3(resid, FR_kerker, size, +1)
    call fft3(resid_cov, FR_wdmetric, size, +1)
    ! deallocate arrays
    deallocate(FR_kerker, FR_wdmetric, STAT=stat)
    if (stat /= 0) &
         call cq_abort("hartree_module/kerker_and_wdmetric: Error &
                        &deallocating FR_kerker, FR_wdmetric: ",  &
                        size, stat)
    call reg_dealloc_mem(area_SC, 2*size, type_dbl)

    return

  end subroutine kerker_and_wdmetric

end module hartree_module
