! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module EarlySC_module
! ------------------------------------------------------------------------------
! Code area 5: self-consistency
! ------------------------------------------------------------------------------

!!****h* Conquest/EarlySC_module
!!  NAME
!!   EarlySC_module
!!  PURPOSE
!!   Collects various routines needed by the early stage search 
!!   for self consistency
!!  USES
!!   datatypes, DiagModule, dimens, DMMin, GenBlas, GenComms, global_module, 
!!   H_matrix_module, logicals, matrix_data, mult_module, numbers
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   28/11/00
!!  MODIFICATION HISTORY
!!   30/01/01 D.R.Bowler
!!    Added brentMin, to implement Brent's method of minimisation
!!   18/05/2001 dave
!!    Stripped down call to get_new_rho in getR2, call to getR2 and
!!    all other calls throughout
!!   06/06/2001 dave
!!    Changes to get_new_rho (use H_matrix_module) and removal of
!!    stop statements. Also added RCS Id and Log tags.
!!   08/06/2001 dave
!!    Converted to use GenComms
!!   17/06/2002 dave
!!    Improved headers slightly, and converted get_new_rho to check
!!    for solution method before generating K (and not generate it if
!!    using diagonalisation)
!!   31/07/2002 dave
!!    Changed matrix_data usage in get_new_rho to use new variable data_M12
!!   2008/02/01 17:47 dave
!!    Changes to write to file not stdout
!!   2012/03/01 L.Tong
!!    Added interfaces for bracketMin, brentMin, get_new_rho, getR2
!!    and reduceLAmbda
!!   2014/09/15 18:30 lat
!!    fixed call start/stop_timer to timer_module (not timer_stdlocks_module !)
!!***
module EarlySCMod

  use datatypes
  use global_module,          only: io_lun, area_SC
  use timer_module,           only: start_timer, stop_timer, cq_timer
  use timer_module,           only: start_backtrace, stop_backtrace
  use timer_stdclocks_module, only: tmr_std_chargescf

  implicit none

  ! Area identification
  integer, parameter, private :: area = 5

  ! RCS ident string for object file id
  character(len=80), save, private :: &
       RCSid = "$Id$"

  !***

contains

  ! -----------------------------------------------------------
  ! Subroutine reduceLambda
  ! -----------------------------------------------------------
  
  !!****f* EarlySC_module/reduceLambda_nospin *
  !!
  !!  NAME 
  !!   reduceLambda_nospin
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Reduces the value of lambda until the residual at that
  !!   point is an acceptable size (with den = lambda den1 + 
  !!   1-lambda den0)
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   28/11/00
  !!  MODIFICATION HISTORY
  !!   18/05/2001 dave
  !!    Stripped call to getR2, indented, added ROBODoc header
  !!   06/06/2001 dave
  !!    Changed stop to cq_abort
  !!   2011/09/13 L.Tong
  !!    Removed obsolete dependence on number_of_bands
  !!   2012/03/01 L.Tong
  !!    Changed name to reduceLambda_nospin, for use with reduceLambda
  !!    interface
  !!   2012/03/18 L.Tong
  !!   - Major rewrite of spin implementation
  !!   - merged subroutines reduceLambda_nospin and reduceLambda_spin
  !!     into one
  !!   - Renamed back to reduceLambda
  !!   - Deleted interface reduceLambda
  !!   - Removed redundant input parameter real(double) mu
  !!  SOURCE
  !!
  subroutine reduceLambda(R0, R1, lambda_1, thresh, MixLin,      &
                          reset_L, fixed_potential, vary_mu,     &
                          n_L_iterations, mx_SCiters, L_tol,     &
                          total_energy, rho0, rho1, resid, size)

    use datatypes
    use numbers
    use GenComms,      only: cq_abort, inode, ionode
    use global_module, only: nspin, spin_factor

    implicit none

    ! Passed variables
    logical      :: MixLin, reset_L, fixed_potential, vary_mu
    integer      :: size
    integer      :: n_L_iterations, mx_SCiters
    real(double) :: R0, R1, lambda_1, thresh
    real(double) :: L_tol, total_energy
    real(double), dimension(:,:) :: rho0, rho1, resid

    ! Local variables
    logical :: reduce
    integer :: nreduc

    reduce = .false.
    if (inode == ionode) &
         write(io_lun, *) 'Welcome to reduceLambda. R0: ', R0
    nreduc = 0
    do while (.not. reduce)
       ! Change lambda and recalculate residual
       nreduc = nreduc + 1
       lambda_1 = half * lambda_1 * (thresh - one) * R0 / (R1 - R0)

       if (inode == ionode) &
            write(io_lun, *) 'reduceLambda: ', nreduc, lambda_1

       R1 = getR2(MixLin, lambda_1, reset_L, fixed_potential, vary_mu, &
                  n_L_iterations, L_tol, total_energy, rho0, rho1,     &
                  resid, size)
       if (R1 <= thresh * R0) then
          reduce = .true.
          if (inode == ionode) write (io_lun, *) 'Done reducing'
       end if
       if (inode == ionode) &
            write (io_lun, *) 'In reduceLambda, R1 is ', R1
       if ((.not. reduce) .and. (nreduc >= mx_SCiters)) then
          call cq_abort('reduceLambda: exceeded reduction iterations',&
                        nreduc, mx_SCiters)
       end if
    end do
    return
  end subroutine reduceLambda
  !!***


  !!****f* EarlySC_module/bracketMin *
  !!
  !!  NAME 
  !!   bracketMin
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Find a triplet of points (a,b,c) such that Rb<Ra and Rb<Rc
  !!   Based on Numerical Recipes in Fortran, 2nd Edition, Section 10.1
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
  !!   18/05/2001 dave
  !!    Reduced subroutine calls
  !!   06/06/2001 dave
  !!    Changed stop to cq_abort
  !!   2011/09/13 L.Tong
  !!    Removed obsolete dependence of number_of_bands
  !!   2011/09/14 L.Tong
  !!    added "::" after "logical" to the line
  !!     logical :: vary_mu, fixed_potential, reset_L
  !!   2012/03/01 L.Tong
  !!    Changed name to bracketMin_nospin, to be used with interface
  !!    bracketMin
  !!   2012/03/18 L.Tong
  !!   - Major rewrite of spin implmentation
  !!   - merged subroutines bracketMin_nospin and bracketMin_spin
  !!   - Renamed to bracketMin
  !!   - deleted interface bracketMin
  !!   - removed redundant input parameter real(double) mu
  !!  SOURCE
  !!
  subroutine bracketMin(moved, lambda_a, lambda_b, lambda_c, Rcross_a, &
                        Rcross_b, Rcross_c, Ra, Rb, Rc, MixLin,        &
                        reset_L, fixed_potential, vary_mu,             &
                        n_L_iterations, mx_SCiters, L_tol,             &
                        total_energy, den0, den1, residb, resid0,      &
                        size)

    use datatypes
    use numbers
    use GenBlas
    use global_module, only: nspin, spin_factor
    use dimens,        only: n_my_grid_points
    use GenComms,      only: cq_abort, gsum, inode, ionode
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    ! Passed variables
    logical :: moved, MixLin ! Flags whether point B has moved
    logical :: vary_mu, fixed_potential, reset_L

    integer :: size
    integer :: mx_SCiters, n_L_iterations

    real(double) :: lambda_a, lambda_b, lambda_c
    real(double) :: Rcross_a, Rcross_b, Rcross_c, Ra, Rb, Rc
    real(double) :: L_tol
    real(double) :: total_energy

    real(double), dimension(:,:) :: den0, den1, residb, resid0

    ! Local variables
    logical      :: bracket
    integer      :: n_brak, spin, stat
    real(double) :: r, cq, u, ulim, Ru, Rcross
    real(double), parameter :: GOLD = 1.618034_double
    real(double), parameter :: GLIMIT = 100.0_double
    real(double), parameter :: TINY = 1.0e-20_double

    real(double), dimension(:,:), allocatable :: residu, residc
    ! I've put residu and residc here as I'm sure that they can share 
    ! space with the n-dimensional residual storage of late-stage

    allocate(residu(size,nspin), residc(size,nspin), STAT=stat)
    if (stat /= 0) call cq_abort("Error allocating residu/c: ", size)
    call reg_alloc_mem(area_SC, 2*nspin*size, type_dbl)

    n_brak = 1
    ! Find an initial third point
    lambda_c = lambda_b + GOLD * (lambda_b - lambda_a)

    Rc = getR2(MixLin, lambda_c, reset_L, fixed_potential, vary_mu, &
               n_L_iterations, L_tol, total_energy, den0, den1,     &
               residc, size)
    
    Rcross_c = zero
    do spin = 1, nspin
       Rcross_c = Rcross_c + spin_factor * &
                  dot(n_my_grid_points, residc(:,spin), 1, resid0(:,spin), 1)
    end do
    ! cross terms
    Rcross_c = Rcross_c + dot(n_my_grid_points, residc(:,1), 1, resid0(:,nspin), 1)
    Rcross_c = Rcross_c + dot(n_my_grid_points, residc(:,nspin), 1, resid0(:,1), 1)
    call gsum (Rcross_c)

    if (inode == ionode) write (io_lun, 208) n_brak, Rc
    
    bracket = .false.
    do while ((.not. bracket) .and. (n_brak < mx_SCiters))
       if (Rc > Rb) then  ! Assume that a>b to start ! Min: (a,b,c)
          bracket = .true.
          if (inode == ionode) write (io_lun, *) 'Accepted bracket'
       else ! Search for a min
          r = (lambda_b - lambda_a) * (Rb - Rc)
          cq = (lambda_b - lambda_c) * (Rb - Ra)
          ! Parabolic extrapolation of u from a,b,c
          u = lambda_b - ((lambda_b - lambda_c) * cq - &
                          (lambda_b - lambda_a) * r) / &
                         (two * sign(max(abs(cq - r), TINY), cq-r))
          ulim = lambda_b + GLIMIT * (lambda_c - lambda_b) ! Limit
          if ((lambda_b - u) * (u - lambda_c) > zero) then  ! b< u < c
             n_brak = n_brak + 1

             Ru = getR2(MixLin, u, reset_L, fixed_potential, vary_mu, &
                        n_L_iterations, L_tol, total_energy, den0,    &
                        den1, residu, size)

             Rcross = zero
             do spin = 1, nspin
                Rcross = Rcross + spin_factor * &
                         dot(n_my_grid_points, residu(:,spin), 1, resid0(:,spin), 1)
             end do
             ! cross terms
             Rcross = Rcross + dot(n_my_grid_points, residu(:,1), 1, &
                                   resid0(:,nspin), 1)
             Rcross = Rcross + dot(n_my_grid_points, residu(:,nspin), 1, &
                                   resid0(:,1), 1)
             call gsum(Rcross)

             if (inode == ionode) write (io_lun, 208) n_brak, Ru

             if (Ru < Rc) then  ! Min between b and c: (b,u,c)
                lambda_a = lambda_b
                Ra = Rb
                Rcross_a = Rcross_b
                lambda_b = u
                moved = .true.
                Rb = Ru
                Rcross_b = Rcross
                do spin = 1, nspin
                   call copy(n_my_grid_points, residu(:,spin), 1, &
                             residb(:,spin), 1)
                end do
                bracket = .true.
             else if (Ru > Rb) then ! Min between a and u: (a,b,u)
                lambda_c = u
                Rc = Ru
                Rcross_c = Rcross
                bracket = .true.
             endif
             if (.not. bracket) then ! Parabolic fit no good -
                !  default 
                u = lambda_c + GOLD * (lambda_c - lambda_b)
                n_brak = n_brak + 1

                Ru = getR2(MixLin, u, reset_L, fixed_potential, &
                           vary_mu, n_L_iterations, L_tol, &
                           total_energy, den0, den1, residu, size)

                if (inode == ionode) write (io_lun, 208) n_brak, Ru

                Rcross = zero
                do spin = 1, nspin
                   Rcross = Rcross + spin_factor * &
                            dot(n_my_grid_points, residu(:,spin), 1, &
                                resid0(:,spin), 1)
                end do
                ! cross terms
                Rcross = Rcross + dot(n_my_grid_points, residu(:,1), &
                                      1, resid0(:,nspin), 1)
                Rcross = Rcross + dot(n_my_grid_points, residu(:,nspin), &
                                      1, resid0(:,1), 1)
                call gsum(Rcross)
             end if

          else if ((lambda_c - u)*(u-ulim) > zero) then ! Fit between c and ulim

             n_brak = n_brak + 1
             
             Ru = getR2(MixLin, u, reset_L, fixed_potential, vary_mu, &
                        n_L_iterations, L_tol, total_energy, den0,    &
                        den1, residu, size)
         
             Rcross = zero
             do spin = 1, nspin
                Rcross = Rcross + spin_factor * &
                         dot(n_my_grid_points, residu(:,spin), 1, resid0(:,spin), 1)
             end do
             ! cross terms
             Rcross = Rcross + dot(n_my_grid_points, residu(:,1), 1, &
                                   resid0(:,nspin), 1)
             Rcross = Rcross + dot(n_my_grid_points, residu(:,nspin), &
                                   1, resid0(:,1), 1)
             call gsum(Rcross)
             
             if (inode == ionode) write (io_lun, 208) n_brak, Ru

             if (Ru < Rc) then  ! Extrapolate u and move all along to (b,c,u)
                lambda_b = lambda_c
                moved = .true.
                lambda_c = u
                u = lambda_c + GOLD * (lambda_c - lambda_b)
                Rb = Rc
                Rc = Ru
                Rcross_b = Rcross_c
                Rcross_c = Rcross
                do spin = 1, nspin
                   residb(1:n_my_grid_points,spin) = residc(1:n_my_grid_points,spin)
                   residc(1:n_my_grid_points,spin) = residu(1:n_my_grid_points,spin)
                end do
                n_brak = n_brak + 1

                Ru = getR2(MixLin, u, reset_L, fixed_potential, &
                           vary_mu, n_L_iterations, L_tol,      &
                           total_energy, den0, den1, residu, size)
                
                Rcross = zero
                do spin = 1, nspin
                   Rcross = Rcross + spin_factor * &
                            dot(n_my_grid_points, residu(:,spin), 1, &
                                resid0(:,spin), 1)
                end do
                ! cross terms
                Rcross = Rcross + dot(n_my_grid_points, residu(:,1), &
                                      1, resid0(:,nspin), 1)
                Rcross = Rcross + dot(n_my_grid_points, residu(:,nspin), &
                                      1, resid0(:,1), 1)
                call gsum(Rcross)

                if (inode == ionode) write (io_lun, 208) n_brak, Ru

             endif

          else if ((u-ulim)*(ulim-lambda_c) >= zero) then ! Parabolic u to ulim

             u = ulim
             n_brak = n_brak + 1

             Ru = getR2(MixLin, u, reset_L, fixed_potential, vary_mu, &
                        n_L_iterations, L_tol, total_energy, den0,    &
                        den1, residu, size)
             
             Rcross = zero
             do spin = 1, nspin
                Rcross = Rcross + spin_factor * &
                         dot(n_my_grid_points, residu(:,spin), 1, resid0(:,spin), 1)
             end do
             ! cross terms
             Rcross = Rcross + dot(n_my_grid_points, residu(:,1), 1, &
                                   resid0(:,nspin), 1)
             Rcross = Rcross + dot(n_my_grid_points, residu(:,nspin), 1, &
                                   resid0(:,1), 1)
             call gsum(Rcross)

             if (inode == ionode) write (io_lun, 208) n_brak, Ru

          else ! Reject parabolic fit and do default magnify

             u = lambda_c + GOLD * (lambda_c - lambda_b)
             n_brak = n_brak + 1
             
             Ru = getR2(MixLin, u, reset_L, fixed_potential, vary_mu, &
                        n_L_iterations, L_tol, total_energy, den0,    &
                        den1, residu, size)
             
             Rcross = zero
             do spin = 1, nspin
                Rcross = Rcross + spin_factor * &
                         dot(n_my_grid_points, residu(:,spin), 1, resid0(:,spin), 1)
             end do
             ! cross terms
             Rcross = Rcross + dot(n_my_grid_points, residu(:,1), 1, &
                                   resid0(:,nspin), 1)
             Rcross = Rcross + dot(n_my_grid_points, residu(:,nspin), &
                                   1, resid0(:,1), 1)
             call gsum(Rcross)


             if (inode == ionode) write (io_lun, 208) n_brak, Ru

          end if

          if (.not. bracket) then ! Update everything to (b,c,u)
             lambda_a = lambda_b
             lambda_b = lambda_c
             lambda_c = u
             moved = .true.
             Ra = Rb
             Rb = Rc
             Rc = Ru
             Rcross_a = Rcross_b
             Rcross_b = Rcross_c
             Rcross_c = Rcross
             do spin = 1, nspin
                residb(1:n_my_grid_points,spin) = residc(1:n_my_grid_points,spin)
                residc(1:n_my_grid_points,spin) = residu(1:n_my_grid_points,spin)
             end do
          end if

       end if ! Rc > Rb
    end do ! while

    if (n_brak >= mx_SCiters) then
       call cq_abort ('bracketMin_spin: Too many attempts', n_brak, mx_SCiters)
    endif

    deallocate(residu, residc, STAT=stat)
    if (stat /= 0) call cq_abort("Error deallocating residu/c")
    call reg_dealloc_mem(area_SC, 2*nspin*size, type_dbl)

208 format('Residual after bracket search ',i5,': ',e15.6)

    return
  end subroutine bracketMin
  !!***


  !!****f* EarlySC_module/searchMin *
  !!
  !!  NAME 
  !!   searchMin
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Based on NR golden search algorithm to reduce the residual - but
  !!   once the residual decreases AT ALL we quit, as this is all we're
  !!   trying to do.
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
  !!   18/05/2001 dave
  !!    Added ROBODoc header, shortened subroutine call to getR2 and
  !!    reduced passed variables
  !!   2011/09/13 L.Tong
  !!    removed obsolete dependence on number_of_bands
  !!   2012/03/18 L.Tong
  !!   - Removed redundant input parameter real(double) mu
  !!   - Added spin implementation
  !!  SOURCE
  !!
  subroutine searchMin(lambda_a,lambda_b,lambda_c,Rcross_a, Rcross_c, &
                       Ra, Rb, Rc, MixLin, reset_L, fixed_potential,  &
                       vary_mu, n_L_iterations, mx_SCiters, L_tol,    &
                       total_energy, den0, den1, residb, size)

    use datatypes
    use numbers
    use global_module, only: nspin, spin_factor
    use dimens,        only: n_my_grid_points
    use GenComms,      only: inode, ionode, cq_abort
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    
    implicit none

    ! Passed variables
    logical :: MixLin, reset_L, fixed_potential, vary_mu

    integer :: size
    integer :: mx_SCiters, n_L_iterations

    real(double) :: Ra, Rb, Rc, lambda_a, lambda_b, lambda_c, &
                    Rcross_a, Rcross_c
    real(double) :: L_tol
    real(double) :: total_energy
    real(double), dimension(:,:) :: den0, den1, residb

    ! Local variables
    logical :: done
    integer :: n_s, spin, stat
    real(double), parameter :: R = 0.61803399_double
    real(double), parameter :: C = 1.0_double - R
    real(double) :: x0, x1, x2, x3
    real(double) :: f0, f1, f2, f3, Rl
    real(double), dimension(:,:), allocatable :: resid

    allocate(resid(size,nspin), STAT=stat)
    if (stat /= 0) call cq_abort("Error allocating resid: ", size, nspin)
    call reg_alloc_mem(area_SC, nspin*size, type_dbl)

    x0 = lambda_a
    x3 = lambda_c
    n_s = 1
    ! Start by subdividing the lambda_a -> lambda_c span in golden ratio
    if(abs(Rc-Rb)>abs(Rb-Ra)) then
       x1 = lambda_b
       x2 = R*lambda_b + C*lambda_c ! or b+C*(c-b)
       f1 = Rb
       f2 = getR2(MixLin, x2, reset_L, fixed_potential, vary_mu,   &
                  n_L_iterations, L_tol, total_energy, den0, den1, &
                  resid, size)
       if (inode == ionode) write (io_lun, 108) n_s, f2
       Rl = f2
    else
       x2 = lambda_b
       x1 = R * lambda_b + C * lambda_a ! or b-C*(b-a)
       f2 = Rb
       f1 = getR2(MixLin, x1, reset_L, fixed_potential, vary_mu,   &
                  n_L_iterations, L_tol, total_energy, den0, den1, &
                  resid, size)
       if (inode == ionode) write (io_lun, 108) n_s, f1
       Rl = f1
    end if
    ! Now loop until done
    done = .false.
    do while ((.not. done) .and. (n_s < mx_SCiters))
       n_s = n_s + 1
       ! Based on which is the lower point, adjust ends and find new middle
       if(f2 < f1) then 
          ! lambda_a = x1
          x0 = x1
          x1 = x2
          x2 = R * x1 + C * x3
          ! Ra = f1
          f0 = f1
          f1 = f2
          f2 = getR2(MixLin, x2, reset_L, fixed_potential, vary_mu,   &
                     n_L_iterations, L_tol, total_energy, den0, den1, &
                     resid, size)
          if (inode == ionode) write (io_lun, 108) n_s, f2
          Rl = f2
          if (Rl < Rb) then
             done = .true.
             ! lambda_a = x1
             lambda_b = x2
             ! Ra = f1
             Rb = f2
             do spin = 1, nspin
                residb(1:n_my_grid_points,spin) = resid(1:n_my_grid_points,spin)
             end do
          endif
       else
          ! lambda_c = x2
          x3 = x2
          x2 = x1
          x1 = R * x2 + C * x0
          Rc = f2
          f3 = f2
          f2 = f1
          f1 = getR2(MixLin, x1, reset_L, fixed_potential, vary_mu,   &
                     n_L_iterations, L_tol, total_energy, den0, den1, &
                     resid, size)
          if (inode == ionode) write (io_lun, 108) n_s, f1
          Rl = f1
          if (Rl < Rb) then
             done = .true.
             ! lambda_c = x2
             lambda_b = x1
             ! Rc = f2
             Rb = f1
             do spin = 1, nspin
                residb(1:n_my_grid_points,spin) = resid(1:n_my_grid_points,spin)
             end do
          endif
       endif
    enddo
    
    if(n_s >= mx_SCiters) then
       if (inode == ionode) &
            write (io_lun, *) &
                  'Attempted to reduce residual via golden search, but failed'
    endif

    deallocate(resid, STAT=stat)
    if (stat /= 0) call cq_abort("Error deallocating resid")
    call reg_dealloc_mem(area_SC, nspin*size, type_dbl)

    return

108 format('Residual after golden section ',i5,': ',e15.6)

  end subroutine searchMin
  !!***


  !!****f* EarlySC_module/brentMin *
  !!
  !!  NAME 
  !!   brentMin
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Based on NR brent algorithm which attempts inverse parabolic interpolation
  !!   but defaults to golden section if the function is uncooperative
  !!   As with searchMin, all we're trying to do is to make the residual decrease
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
  !!   18/05/2001 dave
  !!    Added ROBODoc header, shortened subroutine call to getR2 and
  !!    reduced passed variables
  !!   2011/09/13 L.Tong
  !!    Removed obsolete dependence on number_of_bands
  !!   2012/03/01 L.Tong
  !!    Changed name to brentMin_nospin, to be used with brentMin interface
  !!   2012/03/18 L.Tong
  !!   - Major rewrite of spin implementation
  !!   - merged subroutines brentMin_nospin and brentMin_spin
  !!   - renamed subroutine back to brentMin
  !!   - remved redundant input parameter mu
  !!  SOURCE
  !!
  subroutine brentMin(lambda_a, lambda_b, lambda_c, Rcross_a,      &
                      Rcross_c, Ra, Rb, Rc, MixLin, reset_L,       &
                      fixed_potential, vary_mu, n_L_iterations,    &
                      mx_SCiters, L_tol, total_energy, den0, den1, &
                      residb, size)
    use datatypes
    use numbers
    use global_module, only: nspin, spin_factor
    use dimens, only: n_my_grid_points
    use GenComms, only: inode, ionode, cq_abort
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    
    implicit none

    ! Passed variables
    logical :: MixLin, reset_L, fixed_potential, vary_mu

    integer :: size
    integer :: mx_SCiters, n_L_iterations

    real(double) :: Ra, Rb, Rc, lambda_a, lambda_b, lambda_c, &
                    Rcross_a, Rcross_c
    real(double) :: L_tol
    real(double) :: total_energy
    real(double), dimension(:,:) :: den0, den1, residb

    ! Local variables
    integer      :: n_s, spin, stat
    real(double) :: a, b, d, e, etemp, p, q, r, u, v, w, x, xm
    real(double) :: fu, fv, fw, fx, tol1, tol2
    integer,      parameter :: itmax = 100
    real(double), parameter :: C = 0.3819660_double
    real(double), parameter :: zeps = 1.0e-10_double
    real(double), parameter :: tol = 1.0e-7_double
    real(double), dimension(:,:), allocatable :: resid, residx
    logical :: done

    allocate(resid(size,nspin), residx(size,nspin), STAT=stat)
    if (stat /= 0) call cq_abort("Error allocating resid, residx: ", size, nspin)
    call reg_alloc_mem(area_SC, 2*nspin*size, type_dbl)

    done = .false.
    a = min(lambda_a,lambda_c)
    b = max(lambda_a,lambda_c)
    v = lambda_b
    w = v
    x = v
    e = zero
    fx = Rb
    fv = fx
    fw = fx
    n_s = 0

    do while (.not. done)
       xm = half * (a + b)
       n_s = n_s + 1
       if (inode == ionode) &
            write (io_lun, *) 'Iteration ', n_s,' in brentMin, x, Rx: ', x, fx
       tol1 = tol * abs(x) + zeps
       tol2 = two * tol1
       if (fx < Rb) then ! Not part of the original Brent criterion !
          done = .true.
          lambda_b = x
          Rb = fx
          do spin = 1, nspin
             residb(1:n_my_grid_points,spin) = residx(1:n_my_grid_points,spin)
          end do
          exit ! Leave the do loop
       end if
       if (abs (e) > tol1) then ! Trial parabolic fit
          r = (x - w) * (fx - fv)
          q = (x - v) * (fx - fw)
          p = (x - v) * q - (x - w) * r
          q = two * (q - r)
          if (q > 0) p = -p
          q = abs(q)
          etemp = e
          e = d
          ! Test parabolic fit
          if (abs(p) < abs(half * q * etemp) .or. &
              p > q * (a - x) .or. p < q * (b - x)) then
             d = p / q
             u = x + d
             if (u-a  < tol2 .or. b-u < tol2) d = sign(tol1, xm-x)
          else
             if (x >= xm) then
                e = a - x
             else
                e = b - x
             end if
             d = C * e
          endif
       else ! Golden section
          if (x >= xm) then
             e = a - x
          else
             e = b - x
          end if
          d = C * e
       end if
       ! Check that the step isn't too small
       if (abs(d) >= tol1) then
          u = x + d
       else
          u = x + sign(tol1, d)
       endif
       ! Now get the residual
       fu = getR2(MixLin, u, reset_L, fixed_potential, vary_mu,    &
                  n_L_iterations, L_tol, total_energy, den0, den1, &
                  resid, size)
       if (fu <= fx) then ! Update u to x
          if (u >= x) then
             a = x
          else
             b = x
          end if
          v = w
          fv = fw
          w = x
          fw = fx
          x = u
          fx = fu
          do spin = 1, nspin
             residx(1:n_my_grid_points,spin) = resid(1:n_my_grid_points,spin)
          end do
       else ! Move the end points and continue
          if (u < x) then
             a = u
          else
             b = u
          end if
          if (fu <= fw .or. w == x) then
             v = w
             fv = fw
             w = u
             fw = fu
          else if (fw <= fv .or. v == x .or. v == w) then
             v = u
             fv = fu
          end if
       end if
       if (n_s > itmax) done = .true. ! Force quit
    end do
    if (inode == ionode) &
         write (io_lun, *) 'Finishing brentMin with residual ', Rb

    deallocate(resid, residx, STAT=stat)
    if (stat /= 0) call cq_abort("Error dallocating resid, residx")
    call reg_dealloc_mem(area_SC, 2*nspin*size, type_dbl)

    return
  end subroutine brentMin
  !!*****


  !!****f* EarlySC_module/getR2 *
  !!
  !!  NAME 
  !!   getR2
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Returns the squared norm of the residual for the input formed by going 
  !!   to the point lambda.den1 + (1-lambda).den0, using mixtwo and get_new_rho.
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !!   mixtwo, get_new_rho
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !! 
  !!  MODIFICATION HISTORY
  !!   18/05/2001 dave
  !!    Simplified passed variables and shortened call to get_new_rho
  !!   08/06/2001 dave
  !!    Changed dgsum to gsum
  !!   15:26, 02/05/2005 dave 
  !!    Changed residual in line with SelfConsistency notes
  !!   2011/09/13 L.Tong
  !!    Removed dependence on number_of_bands
  !!   2012/03/01 L.Tong
  !!    Changed name to getR2_nospin, for use with interface getR2
  !!   2012/03/18 L.Tong
  !!   - Major rewrite of spin implementation
  !!   - Merged getR2_nospin and getR2_spin
  !!   - Renamed to getR2
  !!   - Deleted interface getR2
  !!   - Removed redundant input parameter real(double) mu
  !!  SOURCE
  !!
  function getR2(linear, lambda, reset_L, fixed_potential, vary_mu, &
                 n_L_iterations, L_tolerance, total_energy, den0,   &
                 den1, resid, size)
    
    use datatypes
    use numbers
    use GenBlas
    use dimens, only: grid_point_volume
    use GenComms, only: gsum, inode, ionode, cq_abort
    use global_module, only: ne_in_cell, nspin, spin_factor
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none
    
    ! Passed variables
    logical :: linear, reset_L, fixed_potential, vary_mu
    integer :: size
    integer :: n_L_iterations
    real(double) :: lambda
    real(double) :: L_tolerance
    real(double) :: total_energy
    real(double), dimension(:,:) :: den0, den1, resid
    ! result
    real(double) :: getR2
    ! Local variables
    integer      :: spin, stat
    real(double) :: R
    real(double), dimension(:,:), allocatable :: rhoin, rhoout

    allocate(rhoin(size,nspin), rhoout(size,nspin), STAT=stat)
    if (stat /= 0) call cq_abort("Error allocating rhoin/out: ", size, nspin)
    call reg_alloc_mem(area_SC, 2*size*nspin, type_dbl)

    ! Mix the charge densities
    do spin = 1, nspin
       call mixtwo(size, linear, lambda, den0(:,spin), den1(:,spin), &
                   rhoin(:,spin))
    end do

    call get_new_rho(.false., reset_L, fixed_potential, vary_mu, &
                     n_L_iterations, L_tolerance, total_energy,  &
                     rhoin, rhoout, size)
    resid = zero
    R = zero
    do spin = 1, nspin
       call axpy(size, one, rhoout(:,spin), 1, resid(:,spin), 1)
       call axpy(size, -one, rhoin(:,spin), 1, resid(:,spin), 1)
       R = R + spin_factor * dot(size, resid(:,spin), 1, resid(:,spin), 1)
    end do
    ! cross term
    R = R + two * dot(size, resid(:,1), 1, resid(:,nspin), 1)
    call gsum(R)
    getR2 = sqrt(grid_point_volume * R) / ne_in_cell

    deallocate(rhoin, rhoout, STAT=stat)
    if (stat /= 0) call cq_abort("Error deallocating rhoin/out")
    call reg_dealloc_mem(area_SC, 2*size*nspin, type_dbl)

  end function getR2
  !!***


  ! -----------------------------------------------------------
  ! Subroutine get_new_rho
  ! -----------------------------------------------------------
  
  !!****f* EarlySC_module/get_new_rho *
  !!
  !!  NAME 
  !!   get_new_rho
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Gets a new charge density when given an input one; 
  !!   we input rho, find the ground-state density matrix and generate rhoout
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   24/03/99
  !!  MODIFICATION HISTORY
  !!   15/05/2001 dave
  !!    Simplified act_on_vector call and added ROBODoc header
  !!   17/05/2001 dave
  !!    Simplified call to get_H_matrix
  !!   18/05/2001 dave
  !!    Reduced argument list
  !!   30/05/2001 dave
  !!    Included into EarlySC_module
  !!   06/06/2001 dave
  !!    Added use H_matrix_module
  !!   17/06/2002 dave
  !!    Added check for K building depending on solution method
  !!   31/07/2002 dave
  !!    Added data_M12 to use list for matrix_data (follow-on from Pulay
  !!    force in exact diagonalisation)
  !!   08:33, 2003/03/12 dave
  !!    Added call to get_energy to calculate energy
  !!   17:27, 10/05/2005 dave 
  !!    Added (commented out - not as real code) lines to allow checking
  !!    of change of band energy as a way of estimating whether to reset
  !!    L or not in linear scaling work; this is not a fix, but a useful
  !!    idea which I don't want to get lost ! Hence, it's commented out
  !!    below
  !!   10:09, 13/02/2006 drb 
  !!    Removed all explicit references to data_ variables and rewrote
  !!    in terms of new matrix routines
  !!   2006/03/06 05:27 dave
  !!    Tidied calls to get_energy and get_H_matrix
  !!   2011/09/13 L.Tong
  !!    Removed obsolete dependence on number_of_bands
  !!   2011/10/06 13:53 dave
  !!    Added call for cDFT
  !!   2012/03/01 L.Tong
  !!    Changed name to get_new_rho_nospin to use with interface
  !!    get_new_rho
  !!   2012/03/18 L.Tong
  !!   - Rewrote for changed spin implementation, now merges
  !!     subroutines get_new_rho_nospin and get_new_rho_spin
  !!   - Renamed the subroutine back to get_new_rho
  !!   - Removed redundant input parameter real(double) mu
  !!   2015/06/08 lat
  !!    - Added experimental backtrace
  !!   2016/08/01 17:30 nakata
  !!    Introduced atomf instead of sf
  !!   2016/08/08 15:30 nakata
  !!    Renamed supportfns -> atomfns
  !!   2017/02/23 dave
  !!    - Changing location of diagon flag from DiagModule to global and name to flag_diagonalisation
  !!  SOURCE
  !!
  subroutine get_new_rho(record, reset_L, fixed_potential, vary_mu,  &
                         n_CG_L_iterations, tolerance, total_energy, &
                         rhoin, rhoout, size, level)

    use datatypes
    use logicals
    use mult_module,       only: LNV_matrix_multiply
    use DMMin,             only: FindMinDM
    use global_module,     only: iprint_SC, atomf, flag_perform_cDFT, &
                                 nspin, spin_factor, flag_diagonalisation
    use H_matrix_module,   only: get_H_matrix
    !use DiagModule,        only: diagon
    use energy,            only: get_energy, flag_check_DFT
    use functions_on_grid, only: atomfns, allocate_temp_fn_on_grid, &
                                 free_temp_fn_on_grid
    use density_module,    only: get_electronic_density
    use GenComms,          only: inode, ionode
    use cdft_module,       only: cdft_min

    implicit none

    ! Passed Variables
    integer, optional :: level
    logical :: record, reset_L, fixed_potential, vary_mu
    integer :: size, temp_supp_fn
    integer :: n_CG_L_iterations
    real(double) :: tolerance, total_energy
    real(double), dimension(:,:) :: rhoin, rhoout
    type(cq_timer) :: backtrace_timer
    integer        :: backtrace_level

    ! Local variables
    real(double) :: start_BE, new_BE, Ltol
    real(double), dimension(nspin) :: electrons, energy_spin
    
!****lat<$
    if (       present(level) ) backtrace_level = level+1
    if ( .not. present(level) ) backtrace_level = -10
    call start_backtrace(t=backtrace_timer,who='get_new_rho',where=area,&
         level=backtrace_level,echo=.true.)
!****lat>$

    call get_H_matrix(.false., fixed_potential, electrons, rhoin, &
         size, backtrace_level)

    if (flag_perform_cDFT) then
       call cdft_min(reset_L, fixed_potential, vary_mu, &
                     n_CG_L_iterations, tolerance, total_energy)
    else
       Ltol = tolerance
       ! Find minimum density matrix
       call stop_timer(tmr_std_chargescf)
       call FindMinDM(n_CG_L_iterations, vary_mu, Ltol, inode, ionode, &
                      reset_L, record, backtrace_level)
       call start_timer(tmr_std_chargescf)

       ! If we're using O(N), we only have L, and we need K
       ! if diagonalisation, we have K
       if (.not. flag_diagonalisation) then
          call LNV_matrix_multiply(electrons, energy_spin, doK,      &
                                   dontM1,  dontM2, dontM3, dontM4,  &
                                   dontphi, dontE, level=backtrace_level)
       end if
       ! Get total energy
       if (flag_check_DFT) then
          call get_energy(total_energy=total_energy,printDFT=.false., &
               level=backtrace_level)
       else
          call get_energy(total_energy=total_energy,level=backtrace_level)
       endif
    end if ! if (flag_perform_cDFT) then

    ! And get the output density
    temp_supp_fn = allocate_temp_fn_on_grid(atomf)
    call stop_timer(tmr_std_chargescf) ! This routine is always call within area 5
    call get_electronic_density(rhoout, electrons, atomfns, &
                                temp_supp_fn, inode, ionode, size, backtrace_level)
    call start_timer(tmr_std_chargescf)
    call free_temp_fn_on_grid(temp_supp_fn)

    !For DFT energy with charge density constructed by density matrix --
    !TM Nov2007
    if (flag_check_DFT) then
       call get_H_matrix(.false., fixed_potential, electrons, rhoout, size)
       call get_energy(total_energy, flag_check_DFT, backtrace_level)
    endif

!****lat<$
    call stop_backtrace(t=backtrace_timer,who='get_new_rho',echo=.true.)
!****lat>$

    return

  end subroutine get_new_rho
  !*****


  ! -----------------------------------------------------------
  ! Subroutine mixtwo
  ! -----------------------------------------------------------
  
  !!****f* EarlySC_module/mixtwo
  !!
  !!  NAME 
  !!   mixtwo
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Mixes two charge densities (rho0, rho1) according to lambda
  !!   If 0 < lambda < 1, then we get lambda.rho1 + (1-lambda).rho0.
  !!   If lambda < 0 or lambda > 1, then mix in a non-linear manner
  !!   to ensure that the densit is always positive.
  !!   If linear is set to TRUE then linear mixing is always used.
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
  !!   18/05/2001 dave
  !!    Added ROBODoc header and indented
  !!   08/06/2001 dave
  !!    Changed dgsum to gsum
  !!   2012/06/18 L.Tong
  !!   - intialised the sum*s to zero, it used to be non-initialised.
  !!  SOURCE
  !!
  subroutine mixtwo(len, linear, lambda, rho0, rho1, rhomix)

    use datatypes
    use numbers
    use GenBlas
    use GenComms, only: gsum

    implicit none

    ! Passed variables
    integer :: len
    logical :: linear
    real(double) :: lambda
    real(double), dimension(len) :: rho0, rho1, rhomix

    ! Local variables
    integer :: i
    real(double) :: sum0, sum1, summix, factor
    real(double) :: eta, phi, minc, c0, c1, adiffc, sumc
    real(double), parameter :: tol = 0.01_double

    sum0 = zero
    sum1 = zero
    sumc = zero
    summix = zero
    rhomix = zero

    ! Linear mixing
    if(linear.OR.(lambda>zero.AND.lambda<one)) then
       call axpy(len, lambda, rho1, 1, rhomix, 1)
       call axpy(len,(one-lambda), rho0, 1, rhomix, 1)
       ! Non-linear mixing
    else
       do i=1,len  ! Loop over grid points
          c0 = rho0(i)
          c1 = rho1(i)
          sum0 = sum0+c0
          sum1 = sum1+c1
          sumc = c0+c1
          ! Set up useful numbers
          if(c0<c1) then
             minc = c0
             adiffc = c1-c0
             eta = half-lambda
          else
             minc = c1
             adiffc = c0-c1
             eta = lambda-half
          endif
          ! Decide on value of phi
          if(eta<-half) then
             if(minc<=-tol*adiffc*(half+eta))then
                phi = -half*sumc
             else
                phi = -half*sumc + minc*dexp(adiffc*(eta+half)/minc)
             endif
          else if(eta<=half) then
             phi = eta*adiffc
          else
             if(minc<=tol*adiffc*(eta-half)) then
                phi = half*sumc
             else
                phi = half*sumc + minc*dexp(adiffc*(half-eta)/minc)
             endif
          endif ! if(eta<-half)....
          ! Create rhomix
          rhomix(i) = half*sumc+phi
          summix = summix+rhomix(i)
       enddo ! do i=1,len
       ! Normalise
       call gsum(sum0)
       call gsum(sum1)
       call gsum(summix)
       factor = ((one-lambda)*sum0+lambda*sum1)/summix
       call scal(len,factor,rhomix,1)
    endif
    return
  end subroutine mixtwo
  !!***

end module EarlySCMod
