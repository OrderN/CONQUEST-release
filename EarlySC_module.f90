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
!!***
module EarlySCMod

  use datatypes
  use global_module, ONLY: io_lun
  use timer_stdclocks_module, ONLY: start_timer,stop_timer,tmr_std_chargescf

  implicit none

  ! -------------------------------------------------------
  ! RCS ident string for object file id
  ! -------------------------------------------------------
  character(len=80), save, private :: RCSid = "$Id$"

  !!****f* EarlySC_module/bracketMin
  !! PURPOSE
  !!   interface for bracketMin_nospin and bracketMin_spin
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/03/01
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  interface bracketMin
     module procedure bracketMin_nospin
     module procedure bracketMin_spin
  end interface bracketMin
  !!*****

  !!****f* EarlySC_module/brentMin
  !! PURPOSE
  !!   interface for brentMin_nospin and brentMin_spin
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/03/01
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  interface brentMin
     module procedure brentMin_nospin
     module procedure brentMin_spin
  end interface brentMin
  !!*****

  !!****f* EarlySC_module/get_new_rho
  !! PURPOSE
  !!   interface for get_new_rho_nospin and get_new_rho_spin
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/03/01
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  interface get_new_rho
     module procedure get_new_rho_nospin
     module procedure get_new_rho_spin
  end interface get_new_rho
  !!*****

  !!****f* EarlySC_module/getR2
  !! PURPOSE
  !!   interface for getR2_nospin and getR2_spin
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/03/01
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  interface getR2
     module procedure getR2_nospin
     module procedure getR2_spin
  end interface getR2
  !!*****

  !!****f* EarlySC_module/reduceLambda
  !! PURPOSE
  !!   interface for reduceLambda_nospin and reduceLambda_spin
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/03/01
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  interface reduceLambda
     module procedure reduceLambda_nospin
     module procedure reduceLambda_spin
  end interface reduceLambda
  !!*****


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
!!  SOURCE
!!
  subroutine reduceLambda_nospin(R0,R1,lambda_1,thresh,MixLin, reset_L,&
       & fixed_potential, vary_mu, n_L_iterations, mx_SCiters, L_tol,&
       & mu, total_energy, rho0,rho1,resid, size)

    use datatypes
    use numbers
    use GenComms, ONLY: cq_abort, inode, ionode

    implicit none

    ! Passed variables
    logical :: MixLin, reset_L, fixed_potential, vary_mu

    integer :: size
    integer :: n_L_iterations, mx_SCiters

    real(double) :: R0, R1, lambda_1, thresh
    real(double) :: L_tol, mu
    real(double) :: total_energy
    real(double), dimension(size) :: rho0
    real(double), dimension(size) :: rho1
    real(double), dimension(size) :: resid

    ! Local variables
    logical :: reduce
    integer :: nreduc

    reduce = .false.
    if(inode==ionode) write(io_lun,*) 'Welcome to reduceLambda. R0: ',R0
    nreduc = 0
    do while(.NOT.reduce)
       ! Change lambda and recalculate residual
       nreduc = nreduc + 1
       lambda_1 = half*lambda_1*(thresh-one)*R0/(R1-R0)
       if (inode==ionode) write(io_lun,*) 'reduceLambda: ',nreduc,lambda_1
       R1 = getR2( MixLin, lambda_1, reset_L, fixed_potential,&
            & vary_mu, n_L_iterations, L_tol, mu, total_energy, rho0,&
            & rho1,resid, size)
       if(R1<=thresh*R0) then 
          reduce = .true.
          if(inode==ionode) write(io_lun,*) 'Done reducing'
       endif
       if (inode==ionode) write(io_lun,*) 'In reduceLambda, R1 is ',R1
       if((.NOT.reduce).AND.(nreduc>=mx_SCiters)) then
          call cq_abort('reduceLambda: exceeded reduction iterations',&
               nreduc,mx_SCiters)
       endif
    enddo
  end subroutine reduceLambda_nospin
!!***


!!****f* EarlySC_module/reduceLambda_spin
!! PURPOSE
!!   Same as reduceLambda for spin polarised calculations
!!   Based on reduceLambda
!! INPUTS
!! OUTPUT
!! RETURN VALUE
!! AUTHOR
!!   L.Tong
!! CREATION DATE 
!!   2011/09/13
!! MODIFICATION HISTORY
!! SOURCE
!!
  subroutine reduceLambda_spin (R0, R1, lambda_1, thresh, MixLin,&
       & reset_L, fixed_potential, vary_mu, n_L_iterations,&
       & mx_SCiters, L_tol, mu, total_energy, rho0_up, rho0_dn,&
       & rho1_up, rho1_dn, resid_up, resid_dn, size)

    use datatypes
    use numbers
    use GenComms, ONLY: cq_abort, inode, ionode

    implicit none

    ! Passed variables
    logical :: MixLin, reset_L, fixed_potential, vary_mu

    integer :: size
    integer :: n_L_iterations, mx_SCiters

    real(double) :: R0, R1, lambda_1, thresh
    real(double) :: L_tol, mu
    real(double) :: total_energy
    real(double), dimension(size) :: rho0_up, rho0_dn
    real(double), dimension(size) :: rho1_up, rho1_dn
    real(double), dimension(size) :: resid_up, resid_dn

    ! Local variables
    logical :: reduce
    integer :: nreduc

    reduce = .false.
    if (inode == ionode) write(io_lun, *) 'Welcome to reduceLambda.&
         & R0: ',R0
    nreduc = 0
    do while (.not. reduce)
       ! Change lambda and recalculate residual
       nreduc = nreduc + 1
       lambda_1 = half * lambda_1 * (thresh - one) * R0 / (R1-R0)
       if (inode == ionode)&
            & write(io_lun, *) 'reduceLambda: ', nreduc, lambda_1
       R1 = getR2 (MixLin, lambda_1, reset_L, fixed_potential, &
            vary_mu, n_L_iterations, L_tol, mu, total_energy, rho0_up,&
            rho0_dn, rho1_up, rho1_dn, resid_up, resid_dn, size)
       if (R1 <= thresh * R0) then 
          reduce = .true.
          if(inode == ionode) write (io_lun, *) 'Done reducing'
       endif
       if (inode == ionode)&
            & write (io_lun, *) 'In reduceLambda, R1 is ', R1
       if ((.not. reduce) .and. (nreduc >= mx_SCiters)) then
          call cq_abort ('reduceLambda: exceeded reduction iterations',&
               nreduc, mx_SCiters)
       endif
    enddo
  end subroutine reduceLambda_spin
!!*****


!!****f* EarlySC_module/bracketMin_nospin *
!!
!!  NAME 
!!   bracketMin_nospin
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
!!  SOURCE
!!
  subroutine bracketMin_nospin(moved, lambda_a,lambda_b,lambda_c,Rcross_a,&
       & Rcross_b, Rcross_c, Ra, Rb, Rc, MixLin, reset_L,&
       & fixed_potential, vary_mu, n_L_iterations, mx_SCiters, L_tol,&
       & mu, total_energy, den0,den1, residb, resid0, size)

    use datatypes
    use numbers
    use GenBlas
    use dimens, ONLY: n_my_grid_points
    use GenComms, ONLY: cq_abort, gsum, inode, ionode

    implicit none

    ! Passed variables
    logical :: moved, MixLin ! Flags whether point B has moved
    logical :: vary_mu, fixed_potential, reset_L

    integer :: size
    integer :: mx_SCiters, n_L_iterations

    real(double) :: lambda_a, lambda_b, lambda_c
    real(double) :: Rcross_a, Rcross_b, Rcross_c, Ra, Rb, Rc
    real(double) :: L_tol, mu
    real(double) :: total_energy
    real(double), dimension(size) :: den0
    real(double), dimension(size) :: den1
    real(double), dimension(size) :: residb
    real(double), dimension(size) :: resid0

    ! Local variables
    logical :: bracket
    integer :: n_brak
    real(double) :: r, cq, u, ulim, Ru, Rcross
    real(double), parameter :: GOLD = 1.618034_double
    real(double), parameter :: GLIMIT = 100.0_double
    real(double), parameter :: TINY = 1.0e-20_double
    ! Automatic variables
    real(double), dimension(size) :: residu, residc
    ! I've put residu and residc here as I'm sure that they can share 
    ! space with the n-dimensional residual storage of late-stage

    n_brak = 1
    ! Find an initial third point
    lambda_c = lambda_b + GOLD*(lambda_b-lambda_a)
    Rc = getR2( MixLin, lambda_c, reset_L, fixed_potential, vary_mu,&
         & n_L_iterations, L_tol, mu, total_energy, den0,&
         & den1,residc, size)
    Rcross_c = dot(n_my_grid_points, residc, 1, resid0, 1)
    call gsum(Rcross_c)
    if(inode==ionode) write(io_lun,108) n_brak, Rc
108 format('Residual after bracket search ',i5,': ',e15.6)
    bracket = .false.
    do while((.NOT.bracket).AND.(n_brak<mx_SCiters))
       if(Rc>Rb) then  ! Assume that a>b to start ! Min: (a,b,c)
          bracket=.true.
          if(inode==ionode) write(io_lun,*) 'Accepted bracket'
       else ! Search for a min
          r = (lambda_b-lambda_a)*(Rb-Rc)
          cq = (lambda_b - lambda_c)*(Rb-Ra)
          ! Parabolic extrapolation of u from a,b,c
          u = lambda_b - ((lambda_b-lambda_c)*cq-(lambda_b-lambda_a)*r)/ &
               (two*sign(max(abs(cq-r),TINY),cq-r))
          ulim = lambda_b + GLIMIT*(lambda_c-lambda_b) ! Limit
          if((lambda_b-u)*(u-lambda_c)>zero) then  ! b< u < c
             n_brak = n_brak+1
             Ru = getR2( MixLin, u, reset_L, fixed_potential,&
                  & vary_mu, n_L_iterations, L_tol, mu, total_energy,&
                  & den0, den1,residu, size)
             Rcross = dot(n_my_grid_points, residu, 1, resid0, 1)
             call gsum(Rcross)
             if(inode==ionode) write(io_lun,108) n_brak, Ru
             if(Ru<Rc) then  ! Min between b and c: (b,u,c)
                lambda_a = lambda_b
                Ra = Rb
                Rcross_a = Rcross_b
                lambda_b = u
                moved = .true.
                Rb = Ru
                Rcross_b = Rcross
                call copy(n_my_grid_points, residu, 1, residb, 1)
                bracket = .true.
             else if(Ru>Rb) then ! Min between a and u: (a,b,u)
                lambda_c = u
                Rc = Ru
                Rcross_c = Rcross
                bracket = .true.
             endif
             if(.not.bracket) then ! Parabolic fit no good - default 
                u = lambda_c + GOLD*(lambda_c - lambda_b)
                n_brak = n_brak + 1
                Ru = getR2( MixLin, u, reset_L, fixed_potential,&
                     & vary_mu, n_L_iterations, L_tol, mu,&
                     & total_energy, den0, den1,residu, size)
                if(inode==ionode) write(io_lun,108) n_brak, Ru
                Rcross = dot(n_my_grid_points, residu, 1, resid0, 1)
                call gsum(Rcross)
             endif
          else if((lambda_c - u)*(u-ulim)>zero) then ! Fit between c and ulim
             n_brak = n_brak + 1
             Ru = getR2( MixLin, u, reset_L, fixed_potential,&
                  & vary_mu, n_L_iterations, L_tol, mu, total_energy,&
                  & den0, den1,residu, size)
             Rcross = dot(n_my_grid_points, residu, 1, resid0, 1)
             call gsum(Rcross)
             if(inode==ionode) write(io_lun,108) n_brak, Ru
             if(Ru<Rc) then  ! Extrapolate u and move all along to (b,c,u)
                lambda_b = lambda_c
                moved = .true.
                lambda_c = u
                u = lambda_c + GOLD*(lambda_c-lambda_b)
                Rb = Rc
                Rc = Ru
                Rcross_b = Rcross_c
                Rcross_c = Rcross
                residb(1:n_my_grid_points) = residc(1:n_my_grid_points)
                residc(1:n_my_grid_points) = residu(1:n_my_grid_points)
                n_brak = n_brak+1
                Ru = getR2( MixLin, u, reset_L, fixed_potential,&
                     & vary_mu, n_L_iterations, L_tol, mu,&
                     & total_energy, den0, den1,residu, size)
                Rcross = dot(n_my_grid_points, residu, 1, resid0, 1)
                call gsum(Rcross)
                if(inode==ionode) write(io_lun,108) n_brak, Ru
             endif
          else if((u-ulim)*(ulim-lambda_c)>=zero) then ! Parabolic u to ulim
             u = ulim
             n_brak = n_brak +1
             Ru = getR2( MixLin, u, reset_L, fixed_potential,&
                  & vary_mu, n_L_iterations, L_tol, mu, total_energy,&
                  & den0, den1,residu, size)
             Rcross = dot(n_my_grid_points, residu, 1, resid0, 1)
             call gsum(Rcross)
             if(inode==ionode) write(io_lun,108) n_brak, Ru
          else ! Reject parabolic fit and do default magnify
             u = lambda_c + GOLD*(lambda_c-lambda_b)
             n_brak = n_brak + 1
             Ru = getR2( MixLin, u, reset_L, fixed_potential,&
                  & vary_mu, n_L_iterations, L_tol, mu, total_energy,&
                  & den0, den1,residu, size)
             Rcross = dot(n_my_grid_points, residu, 1, resid0, 1)
             call gsum(Rcross)
             if(inode==ionode) write(io_lun,108) n_brak, Ru
          endif
          if(.NOT.bracket) then ! Update everything to (b,c,u)
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
             residb(1:n_my_grid_points) = residc(1:n_my_grid_points)
             residc(1:n_my_grid_points) = residu(1:n_my_grid_points)
          endif
       endif
    enddo
    if(n_brak>=mx_SCiters) then
       call cq_abort('bracketMin: Too many attempts',n_brak,mx_SCiters)
    endif
  end subroutine bracketMin_nospin
!!***


!!****f* EarlySC_module/bracketMin_spin
!! PURPOSE
!!   Same purpose as for bracketMin, but for spin polarisation
!!   calculations. Based on bracketMin
!! INPUTS
!! OUTPUT
!! RETURN VALUE
!! AUTHOR
!!   L.Tong
!! CREATION DATE 
!!   2011/09/14
!! MODIFICATION HISTORY
!! SOURCE
!!
  subroutine bracketMin_spin (moved, lambda_a, lambda_b, lambda_c,&
       & Rcross_a, Rcross_b, Rcross_c, Ra, Rb, Rc, MixLin, reset_L,&
       & fixed_potential, vary_mu, n_L_iterations, mx_SCiters, L_tol,&
       & mu, total_energy, den0_up, den0_dn, den1_up, den1_dn,&
       & residb_up, residb_dn, resid0_up, resid0_dn, size)

    use datatypes
    use numbers
    use GenBlas
    use dimens, ONLY: n_my_grid_points
    use GenComms, ONLY: cq_abort, gsum, inode, ionode

    implicit none

    ! Passed variables
    logical :: moved, MixLin ! Flags whether point B has moved
    logical :: vary_mu, fixed_potential, reset_L

    integer :: size
    integer :: mx_SCiters, n_L_iterations

    real(double) :: lambda_a, lambda_b, lambda_c
    real(double) :: Rcross_a, Rcross_b, Rcross_c, Ra, Rb, Rc
    real(double) :: L_tol, mu
    real(double) :: total_energy
    real(double), dimension(size) :: den0_up, den0_dn
    real(double), dimension(size) :: den1_up, den1_dn
    real(double), dimension(size) :: residb_up, residb_dn
    real(double), dimension(size) :: resid0_up, resid0_dn

    ! Local variables
    logical :: bracket
    integer :: n_brak
    real(double) :: r, cq, u, ulim, Ru, Rcross
    real(double), parameter :: GOLD = 1.618034_double
    real(double), parameter :: GLIMIT = 100.0_double
    real(double), parameter :: TINY = 1.0e-20_double
    ! Automatic variables
    real(double), dimension(size) :: residu_up, residc_up
    real(double), dimension(size) :: residu_dn, residc_dn
    ! I've put residu and residc here as I'm sure that they can share 
    ! space with the n-dimensional residual storage of late-stage

    n_brak = 1
    ! Find an initial third point
    lambda_c = lambda_b + GOLD*(lambda_b-lambda_a)
    Rc = getR2 (MixLin, lambda_c, reset_L, fixed_potential, &
         vary_mu, n_L_iterations, L_tol, mu, total_energy, den0_up, &
         den0_dn, den1_up, den1_dn, residc_up, residc_dn, size)
    Rcross_c = dot (n_my_grid_points, residc_up, 1, resid0_up, 1)
    Rcross_c = Rcross_c + dot (n_my_grid_points, residc_dn, 1, resid0_dn, 1)
    Rcross_c = Rcross_c + dot (n_my_grid_points, residc_up, 1, resid0_dn, 1)
    Rcross_c = Rcross_c + dot (n_my_grid_points, residc_dn, 1, resid0_up, 1)
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
          u = lambda_b - ((lambda_b - lambda_c) * cq - (lambda_b -&
               & lambda_a) * r) / (two*sign(max(abs(cq-r),TINY),cq-r))
          ulim = lambda_b + GLIMIT * (lambda_c - lambda_b) ! Limit
          if ((lambda_b - u) * (u - lambda_c) > zero) then  ! b< u < c
             n_brak = n_brak + 1
             Ru = getR2 (MixLin, u, reset_L, fixed_potential,&
                  & vary_mu, n_L_iterations, L_tol, mu, total_energy,&
                  & den0_up, den0_dn, den1_up, den1_dn, residu_up,&
                  & residu_dn, size)
             Rcross = dot (n_my_grid_points, residu_up, 1, resid0_up, 1)
             Rcross = Rcross + dot (n_my_grid_points, residu_dn, 1, resid0_dn, 1)
             Rcross = Rcross + dot (n_my_grid_points, residu_up, 1, resid0_dn, 1)
             Rcross = Rcross + dot (n_my_grid_points, residu_dn, 1,&
                  & resid0_up, 1)
             call gsum (Rcross)
             if (inode == ionode) write (io_lun, 208) n_brak, Ru
             if (Ru < Rc) then  ! Min between b and c: (b,u,c)
                lambda_a = lambda_b
                Ra = Rb
                Rcross_a = Rcross_b
                lambda_b = u
                moved = .true.
                Rb = Ru
                Rcross_b = Rcross
                call copy (n_my_grid_points, residu_up, 1, residb_up,&
                     & 1)
                call copy (n_my_grid_points, residu_dn, 1, residb_dn,&
                     & 1)
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
                Ru = getR2 (MixLin, u, reset_L, fixed_potential,&
                     & vary_mu, n_L_iterations, L_tol, mu,&
                     & total_energy, den0_up, den0_dn, den1_up,&
                     & den1_dn, residu_up, residu_dn,  size)
                if (inode==ionode) write (io_lun, 208) n_brak, Ru
                Rcross = dot (n_my_grid_points, residu_up, 1,&
                     & resid0_up, 1)
                Rcross = Rcross + dot (n_my_grid_points, residu_dn,&
                     & 1, resid0_dn, 1)
                Rcross = Rcross + dot (n_my_grid_points, residu_up,&
                     & 1, resid0_dn, 1)
                Rcross = Rcross + dot (n_my_grid_points, residu_dn,&
                     & 1, resid0_up, 1)
                call gsum (Rcross)
             endif
          else if ((lambda_c - u)*(u-ulim) > zero) then ! Fit between c and ulim
             n_brak = n_brak + 1
             Ru = getR2 (MixLin, u, reset_L, fixed_potential,&
                  & vary_mu, n_L_iterations, L_tol, mu, total_energy,&
                  & den0_up, den0_dn, den1_up, den1_dn, residu_up,&
                  & residu_dn, size)
             Rcross = dot (n_my_grid_points, residu_up, 1, resid0_up, 1)
             Rcross = Rcross + dot (n_my_grid_points, residu_dn, 1, resid0_dn, 1)
             Rcross = Rcross + dot (n_my_grid_points, residu_up, 1, resid0_dn, 1)
             Rcross = Rcross + dot (n_my_grid_points, residu_dn, 1, resid0_up, 1)
             call gsum (Rcross)
             if (inode == ionode) write (io_lun, 208) n_brak, Ru
             if (Ru < Rc) then  ! Extrapolate u and move all along to (b,c,u)
                lambda_b = lambda_c
                moved = .true.
                lambda_c = u
                u = lambda_c + GOLD*(lambda_c-lambda_b)
                Rb = Rc
                Rc = Ru
                Rcross_b = Rcross_c
                Rcross_c = Rcross
                residb_up(1:n_my_grid_points) = residc_up(1:n_my_grid_points)
                residb_dn(1:n_my_grid_points) = residc_dn(1:n_my_grid_points)
                residc_up(1:n_my_grid_points) = residu_up(1:n_my_grid_points)
                residc_dn(1:n_my_grid_points) = residu_dn(1:n_my_grid_points)
                n_brak = n_brak + 1
                Ru = getR2 (MixLin, u, reset_L, fixed_potential,&
                     & vary_mu, n_L_iterations, L_tol, mu,&
                     & total_energy, den0_up, den0_dn, den1_up,&
                     & den1_dn, residu_up, residu_dn, size)
                Rcross = dot (n_my_grid_points, residu_up, 1, resid0_up, 1)
                Rcross = Rcross + dot (n_my_grid_points, residu_dn,&
                     & 1, resid0_dn, 1)
                Rcross = Rcross + dot (n_my_grid_points, residu_up,&
                     & 1, resid0_dn, 1)
                Rcross = Rcross + dot (n_my_grid_points, residu_dn,&
                     & 1, resid0_up, 1)
                call gsum (Rcross)
                if (inode == ionode) write (io_lun, 208) n_brak, Ru
             endif
          else if ((u-ulim)*(ulim-lambda_c) >= zero) then ! Parabolic u to ulim
             u = ulim
             n_brak = n_brak + 1
             Ru = getR2 (MixLin, u, reset_L, fixed_potential,&
                  & vary_mu, n_L_iterations, L_tol, mu, total_energy,&
                  & den0_up, den0_dn, den1_up, den1_dn, residu_up,&
                  & residu_dn, size)
             Rcross = dot (n_my_grid_points, residu_up, 1, resid0_up, 1)
             Rcross = Rcross + dot (n_my_grid_points, residu_dn, 1, resid0_dn, 1)
             Rcross = Rcross + dot (n_my_grid_points, residu_up, 1, resid0_dn, 1)
             Rcross = Rcross + dot (n_my_grid_points, residu_dn, 1, resid0_up, 1)
             call gsum (Rcross)
             if (inode == ionode) write (io_lun, 208) n_brak, Ru
          else ! Reject parabolic fit and do default magnify
             u = lambda_c + GOLD * (lambda_c-lambda_b)
             n_brak = n_brak + 1
             Ru = getR2 (MixLin, u, reset_L, fixed_potential,&
                  & vary_mu, n_L_iterations, L_tol, mu, total_energy,&
                  & den0_up, den0_dn, den1_up, den1_dn, residu_up,&
                  & residu_dn, size)
             Rcross = dot (n_my_grid_points, residu_up, 1, resid0_up, 1)
             Rcross = Rcross + dot (n_my_grid_points, residu_dn, 1, resid0_dn, 1)
             Rcross = Rcross + dot (n_my_grid_points, residu_up, 1, resid0_dn, 1)
             Rcross = Rcross + dot (n_my_grid_points, residu_dn, 1, resid0_up, 1)
             call gsum (Rcross)
             if (inode == ionode) write (io_lun, 208) n_brak, Ru
          endif
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
             residb_up(1:n_my_grid_points) = residc_up(1:n_my_grid_points)
             residb_dn(1:n_my_grid_points) = residc_dn(1:n_my_grid_points)
             residc_up(1:n_my_grid_points) = residu_up(1:n_my_grid_points)
             residc_dn(1:n_my_grid_points) = residu_dn(1:n_my_grid_points)
          endif
       endif
    enddo

    if (n_brak >= mx_SCiters) then
       call cq_abort ('bracketMin_spin: Too many attempts', n_brak, mx_SCiters)
    endif

208 format('Residual after bracket search ',i5,': ',e15.6)

    return
 end subroutine bracketMin_spin
!!*****


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
!!  SOURCE
!!
  subroutine searchMin(lambda_a,lambda_b,lambda_c,Rcross_a, Rcross_c,&
       & Ra, Rb, Rc, MixLin, reset_L, fixed_potential, vary_mu,&
       & n_L_iterations, mx_SCiters, L_tol, mu, total_energy,&
       & den0,den1, residb, size)

    use datatypes
    use numbers
    use dimens, ONLY: n_my_grid_points
    use GenComms, ONLY: inode, ionode
    
    implicit none

    ! Passed variables
    logical :: MixLin, reset_L, fixed_potential, vary_mu

    integer :: size
    integer :: mx_SCiters, n_L_iterations

    real(double) :: Ra, Rb, Rc, lambda_a, lambda_b, lambda_c, &
         Rcross_a, Rcross_c
    real(double) :: L_tol, mu
    real(double) :: total_energy
    real(double), dimension(size) :: den0
    real(double), dimension(size) :: den1
    real(double), dimension(size) :: residb

    ! Local variables
    logical :: done
    integer :: n_s
    real(double), parameter :: R = 0.61803399_double
    real(double), parameter :: C = 1.0_double - R
    real(double) :: x0, x1, x2, x3
    real(double) :: f0, f1, f2, f3, Rl
    real(double), dimension(size) :: resid, rhoin, rhoout

    x0 = lambda_a
    x3 = lambda_c
    n_s = 1
    ! Start by subdividing the lambda_a -> lambda_c span in golden ratio
    if(abs(Rc-Rb)>abs(Rb-Ra)) then
       x1 = lambda_b
       x2 = R*lambda_b + C*lambda_c ! or b+C*(c-b)
       f1 = Rb
       f2 = getR2( MixLin, x2, reset_L, fixed_potential, vary_mu,&
            & n_L_iterations, L_tol, mu, total_energy, den0,&
            & den1,resid, size)
       if(inode==ionode) write(io_lun,108) n_s, f2
       Rl = f2
    else
       x2 = lambda_b
       x1 = R*lambda_b + C*lambda_a ! or b-C*(b-a)
       f2 = Rb
       f1 = getR2( MixLin, x1, reset_L, fixed_potential, vary_mu,&
            & n_L_iterations, L_tol, mu, total_energy, den0,&
            & den1,resid, size)
       if(inode==ionode) write(io_lun,108) n_s, f1
       Rl = f1
    endif
    ! Now loop until done
    done = .false.
    do while((.NOT.done).AND.(n_s<mx_SCiters))
       n_s = n_s + 1
       ! Based on which is the lower point, adjust ends and find new middle
       if(f2 < f1) then 
          !         lambda_a = x1
          x0 = x1
          x1 = x2
          x2 = R*x1+C*x3
          !         Ra = f1
          f0 = f1
          f1 = f2
          f2 = getR2( MixLin, x2, reset_L, fixed_potential, vary_mu,&
               & n_L_iterations, L_tol, mu, total_energy, den0,&
               & den1,resid, size)
          if(inode==ionode) write(io_lun,108) n_s, f2
          Rl = f2
          if(Rl<Rb) then
             done = .true.
             !            lambda_a = x1
             lambda_b = x2
             !            Ra = f1
             Rb = f2
             residb(1:n_my_grid_points) = resid(1:n_my_grid_points)
          endif
       else
          !         lambda_c = x2
          x3 = x2
          x2 = x1
          x1 = R*x2 + C*x0
          Rc = f2
          f3 = f2
          f2 = f1
          f1 = getR2( MixLin, x1, reset_L, fixed_potential, vary_mu,&
               & n_L_iterations, L_tol, mu, total_energy, den0,&
               & den1,resid, size)
          if(inode==ionode) write(io_lun,108) n_s, f1
          Rl = f1
          if(Rl<Rb) then
             done = .true.
             !            lambda_c = x2
             lambda_b = x1
             !            Rc = f2
             Rb = f1
             residb(1:n_my_grid_points) = resid(1:n_my_grid_points)
          endif
       endif
    enddo
108 format('Residual after golden section ',i5,': ',e15.6)
    if(n_s>=mx_SCiters) then
       if(inode==ionode) write(io_lun,*) 'Attempted to reduce residual via golden search, but failed'
    endif
  end subroutine searchMin
!!***

!!****f* EarlySC_module/brentMin_nospin *
!!
!!  NAME 
!!   brentMin_nospin
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
!!  SOURCE
!!
  subroutine brentMin_nospin(lambda_a,lambda_b,lambda_c,Rcross_a, Rcross_c,&
       & Ra, Rb, Rc, MixLin, reset_L, fixed_potential, vary_mu,&
       & n_L_iterations, mx_SCiters, L_tol, mu, total_energy,&
       & den0,den1, residb, size)

    use datatypes
    use numbers
    use dimens, ONLY: n_my_grid_points
    use GenComms, ONLY: inode, ionode
    
    implicit none

    ! Passed variables
    logical :: MixLin, reset_L, fixed_potential, vary_mu

    integer :: size
    integer :: mx_SCiters, n_L_iterations

    real(double) :: Ra, Rb, Rc, lambda_a, lambda_b, lambda_c,&
         & Rcross_a, Rcross_c
    real(double) :: L_tol, mu
    real(double) :: total_energy
    real(double), dimension(size) :: den0
    real(double), dimension(size) :: den1
    real(double), dimension(size) :: residb

    ! Local variables
    integer :: n_s
    integer, parameter :: itmax = 100
    real(double), parameter :: C = 0.3819660_double
    real(double), parameter :: zeps = 1.0e-10_double
    real(double), parameter :: tol = 1.0e-7_double
    real(double) :: a,b,d,e,etemp,p,q,r,u,v,w,x,xm
    real(double) :: fu,fv,fw,fx,tol1,tol2
    real(double), dimension(size) :: resid, residx
    logical :: done

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
    do while(.NOT.done)
       xm = half*(a+b)
       n_s = n_s+1
       if(inode==ionode) write(io_lun,*) 'Iteration ',n_s,' in&
            & brentMin, x,Rx: ',x,fx
       tol1 = tol*abs(x)+zeps
       tol2 = two*tol1
       if(fx<Rb) then ! Not part of the original Brent criterion !
          done = .true.
          lambda_b = x
          Rb = fx
          residb(1:n_my_grid_points) = residx(1:n_my_grid_points)
          exit ! Leave the do loop
       endif
       if(abs(e)>tol1) then ! Trial parabolic fit
          r = (x-w)*(fx-fv)
          q = (x-v)*(fx-fw)
          p = (x-v)*q-(x-w)*r
          q = two*(q-r)
          if(q>0) p=-p
          q = abs(q)
          etemp = e
          e = d
          ! Test parabolic fit
          if(abs(p)<abs(half*q*etemp).OR.p>q*(a-x).OR.p<q*(b-x)) then
             d=p/q
             u=x+d
             if(u-a<tol2.OR.b-u<tol2) d=sign(tol1,xm-x)
          else
             if(x>=xm) then
                e=a-x
             else
                e=b-x
             endif
             d = C*e
          endif
       else ! Golden section
          if(x>=xm) then
             e=a-x
          else
             e=b-x
          endif
          d = C*e
       endif
       ! Check that the step isn't too small
       if(abs(d)>=tol1) then
          u=x+d
       else
          u=x+sign(tol1,d)
       endif
       ! Now get the residual
       fu = getR2( MixLin, u, reset_L, fixed_potential, vary_mu,&
            & n_L_iterations, L_tol, mu, total_energy, den0,&
            & den1,resid, size)
       if(fu<=fx) then ! Update u to x
          if(u>=x) then
             a = x
          else
             b = x
          endif
          v=w
          fv = fw
          w=x
          fw = fx
          x=u
          fx = fu
          residx(1:n_my_grid_points) = resid(1:n_my_grid_points)
       else ! Move the end points and continue
          if(u<x) then
             a = u
          else
             b = u
          endif
          if(fu<=fw.OR.w==x) then
             v=w
             fv = fw
             w=u
             fw = fu
          else if(fw<=fv.OR.v==x.OR.v==w) then
             v=u
             fv = fu
          endif
       endif
       if(n_s>itmax) done = .true. ! Force quit
    enddo
    if(inode==ionode) write(io_lun,*) 'Finishing brentMin with residual ',Rb
    return
  end subroutine brentMin_nospin
!!***


!!****f* EarlySC_module/brentMin_spin
!! PURPOSE
!!   Same as brentMin, but for spin polarised calculations.
!!   Based on brentMin.
!! INPUTS
!! OUTPUT
!! RETURN VALUE
!! AUTHOR
!!   L.Tong
!! CREATION DATE 
!!   2011/09/14
!! MODIFICATION HISTORY
!! SOURCE
!!
  subroutine brentMin_spin (lambda_a, lambda_b, lambda_c, Rcross_a,&
       & Rcross_c, Ra, Rb, Rc, MixLin, reset_L, fixed_potential,&
       & vary_mu, n_L_iterations, mx_SCiters, L_tol, mu,&
       & total_energy, den0_up, den0_dn, den1_up, den1_dn, residb_up,&
       & resid_dn, size)

    use datatypes
    use numbers
    use dimens, ONLY: n_my_grid_points
    use GenComms, ONLY: inode, ionode
    
    implicit none

    ! Passed variables
    logical :: MixLin, reset_L, fixed_potential, vary_mu

    integer :: size
    integer :: mx_SCiters, n_L_iterations

    real(double) :: Ra, Rb, Rc, lambda_a, lambda_b, lambda_c,&
         & Rcross_a, Rcross_c
    real(double) :: L_tol, mu
    real(double) :: total_energy
    real(double), dimension(size) :: den0_up, den0_dn
    real(double), dimension(size) :: den1_up, den1_dn
    real(double), dimension(size) :: residb_up, residb_dn

    ! Local variables
    integer :: n_s
    integer, parameter :: itmax = 100
    real(double), parameter :: C = 0.3819660_double
    real(double), parameter :: zeps = 1.0e-10_double
    real(double), parameter :: tol = 1.0e-7_double
    real(double) :: a,b,d,e,etemp,p,q,r,u,v,w,x,xm
    real(double) :: fu,fv,fw,fx,tol1,tol2
    real(double), dimension(size) :: resid_up, residx_up
    real(double), dimension(size) :: resid_dn, residx_dn
    logical :: done

    done = .false.
    a = min (lambda_a,lambda_c)
    b = max (lambda_a,lambda_c)
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
       if (inode == ionode) write (io_lun, *) 'Iteration ', n_s,' in&
            & brentMin, x,Rx: ', x, fx
       tol1 = tol * abs (x) + zeps
       tol2 = two * tol1
       if (fx < Rb) then ! Not part of the original Brent criterion !
          done = .true.
          lambda_b = x
          Rb = fx
          residb_up(1:n_my_grid_points) = residx_up(1:n_my_grid_points)
          residb_dn(1:n_my_grid_points) = residx_dn(1:n_my_grid_points)
          exit ! Leave the do loop
       end if
       if (abs (e) > tol1) then ! Trial parabolic fit
          r = (x - w) * (fx - fv)
          q = (x - v) * (fx - fw)
          p = (x - v) * q - (x - w) * r
          q = two * (q - r)
          if (q > 0) p = -p
          q = abs (q)
          etemp = e
          e = d
          ! Test parabolic fit
          if (abs (p) < abs (half * q * etemp) .or. p > q * (a - x)&
               & .or. p < q * (b - x)) then
             d = p / q
             u = x + d
             if (u-a  < tol2 .or. b-u < tol2) d = sign (tol1, xm-x)
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
          endif
          d = C * e
       endif
       ! Check that the step isn't too small
       if (abs (d) >= tol1) then
          u = x + d
       else
          u = x + sign (tol1, d)
       endif
       ! Now get the residual
       fu = getR2 (MixLin, u, reset_L, fixed_potential, vary_mu, &
            n_L_iterations, L_tol, mu, total_energy, den0_up, den0_dn, &
            den1_up, den1_dn, resid_up, resid_dn, size)
       if (fu <= fx) then ! Update u to x
          if (u >= x) then
             a = x
          else
             b = x
          endif
          v = w
          fv = fw
          w = x
          fw = fx
          x = u
          fx = fu
          residx_up(1:n_my_grid_points) = resid_up(1:n_my_grid_points)
          residx_dn(1:n_my_grid_points) = resid_dn(1:n_my_grid_points)
       else ! Move the end points and continue
          if (u < x) then
             a = u
          else
             b = u
          endif
          if (fu <= fw .or. w == x) then
             v = w
             fv = fw
             w = u
             fw = fu
          else if (fw <= fv .or. v == x .or. v == w) then
             v = u
             fv = fu
          endif
       endif
       if (n_s > itmax) done = .true. ! Force quit
    enddo
    if (inode == ionode)&
         & write (io_lun, *) 'Finishing brentMin with residual ', Rb
    return
  end subroutine brentMin_spin
!!*****


!!****f* EarlySC_module/getR2_nospin *
!!
!!  NAME 
!!   getR2_nospin
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
!!  SOURCE
!!
  function getR2_nospin(linear, lambda, reset_L, fixed_potential, &
       vary_mu, n_L_iterations, L_tolerance, mu, total_energy, den0, &
       den1,resid, size)

    use datatypes
    use numbers
    use GenBlas
    use dimens, ONLY: grid_point_volume
    use GenComms, ONLY: gsum, inode, ionode
    use global_module, ONLY: ne_in_cell

    implicit none

    ! Passed variables
    logical :: linear, reset_L, fixed_potential, vary_mu
    ! result
    real(double) :: getR2_nospin

    integer :: size
    integer :: n_L_iterations

    real(double) :: lambda
    real(double) :: L_tolerance, mu
    real(double) :: total_energy
    real(double), dimension(size) :: den0
    real(double), dimension(size) :: den1
    real(double), dimension(size) :: resid

    ! Local variables
    real(double) :: R
    ! Automatic arrays
    real(double), dimension(size) :: rhoin, rhoout

    ! Mix the charge densities
    call mixtwo(size, linear, lambda, den0, den1, rhoin)
    call get_new_rho( .false., reset_L, fixed_potential, vary_mu, &
         n_L_iterations, L_tolerance, mu, total_energy, rhoin, rhoout,&
         size)
    resid = zero
    call axpy(size, one, rhoout, 1, resid, 1)
    call axpy(size, -one, rhoin, 1, resid, 1)
    R = dot(size,resid,1,resid,1)
    call gsum(R)
    getR2_nospin = sqrt(grid_point_volume*R)/ne_in_cell
  end function getR2_nospin
!!***


!!****f* EarlySC_module/getR2_spin
!! PURPOSE
!!   Same as getR2, but for spin polarised calculations
!! INPUTS
!! OUTPUT
!! RETURN VALUE
!! AUTHOR
!!   L.Tong
!! CREATION DATE 
!!   2011/09/13
!! MODIFICATION HISTORY
!! SOURCE
!!
  function getR2_spin (linear, lambda, reset_L, fixed_potential, &
       vary_mu, n_L_iterations, L_tolerance, mu, total_energy, &
       den0_up, den0_dn, den1_up, den1_dn, resid_up, resid_dn, size)
    
    use datatypes
    use numbers
    use GenBlas
    use dimens, ONLY: grid_point_volume
    use GenComms, ONLY: gsum, inode, ionode
    use global_module, ONLY: ne_in_cell

    implicit none
    
    ! Passed variables
    logical :: linear, reset_L, fixed_potential, vary_mu
    integer :: size
    integer :: n_L_iterations
    ! result
    real(double) :: getR2_spin
  
    real(double) :: lambda
    real(double) :: L_tolerance, mu
    real(double) :: total_energy
    real(double), dimension(size) :: den0_up, den0_dn
    real(double), dimension(size) :: den1_up, den1_dn
    real(double), dimension(size) :: resid_up, resid_dn

    ! Local variables
    real(double) :: R
    ! Automatic arrays
    real(double), dimension(size) :: rhoin_up, rhoin_dn, rhoout_up,&
         & rhoout_dn

    ! Mix the charge densities
    call mixtwo (size, linear, lambda, den0_up, den1_up, rhoin_up)
    call mixtwo (size, linear, lambda, den0_dn, den1_dn, rhoin_dn)
    call get_new_rho (.false., reset_L, fixed_potential, vary_mu,&
         n_L_iterations, L_tolerance, mu, total_energy, rhoin_up, &
         rhoin_dn, rhoout_up, rhoout_dn, size)
    resid_up = zero
    resid_dn = zero
    call axpy (size, one, rhoout_up, 1, resid_up, 1)
    call axpy (size, -one, rhoin_dn, 1, resid_dn, 1)
    R = dot (size, resid_up, 1, resid_up, 1)
    R = R + dot (size, resid_dn, 1, resid_dn, 1)
    R = R + two * dot (size, resid_up, 1, resid_dn, 1)
    call gsum(R)
    getR2_spin = sqrt (grid_point_volume * R) / ne_in_cell

  end function getR2_spin
!!*****


! -----------------------------------------------------------
! Subroutine get_new_rho
! -----------------------------------------------------------

!!****f* EarlySC_module/get_new_rho_nospin *
!!
!!  NAME 
!!   get_new_rho_nospin
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
!!    Added data_M12 to use list for matrix_data (follow-on from Pulay force in exact diagonalisation)
!!   08:33, 2003/03/12 dave
!!    Added call to get_energy to calculate energy
!!   17:27, 10/05/2005 dave 
!!    Added (commented out - not as real code) lines to allow checking of change of band energy as a way
!!    of estimating whether to reset L or not in linear scaling work; this is not a fix, but a useful idea
!!    which I don't want to get lost ! Hence, it's commented out below
!!   10:09, 13/02/2006 drb 
!!    Removed all explicit references to data_ variables and rewrote in terms of new 
!!    matrix routines
!!   2006/03/06 05:27 dave
!!    Tidied calls to get_energy and get_H_matrix
!!   2011/09/13 L.Tong
!!    Removed obsolete dependence on number_of_bands
!!   2011/10/06 13:53 dave
!!    Added call for cDFT
!!   2012/03/01 L.Tong
!!    Changed name to get_new_rho_nospin to use with interface get_new_rho
!!  SOURCE
!!
subroutine get_new_rho_nospin( record, reset_L, fixed_potential, vary_mu,&
     & n_CG_L_iterations, tolerance, mu, total_energy, rhoin, rhoout,&
     & size)
  
  use datatypes
  use logicals
  use mult_module, ONLY: LNV_matrix_multiply
  use DMMin, ONLY: FindMinDM
  use global_module, ONLY: iprint_SC, sf, flag_perform_cDFT
  use H_matrix_module, ONLY: get_H_matrix
  use DiagModule, ONLY: diagon
  use energy, ONLY: get_energy, flag_check_DFT
  use functions_on_grid, ONLY: supportfns, allocate_temp_fn_on_grid, free_temp_fn_on_grid
  use density_module, ONLY: get_electronic_density
  use GenComms, ONLY: inode, ionode
  use cdft_module, ONLY: cdft_min

  implicit none
  
  ! Passed Variables
  logical :: record, reset_L, fixed_potential, vary_mu
  
  integer :: size, temp_supp_fn
  integer :: n_CG_L_iterations
  
  real(double) :: tolerance, mu, total_energy
  real(double), dimension(size) :: rhoin
  real(double), dimension(size) :: rhoout
  
  ! Local variables
  real(double) :: electrons, total_energy_1, start_BE, new_BE, Ltol
  !For DFT energy with charge density constructed by density matrix   --  TM  Nov2007

  ! ** Useful, but not rigorous **
  !start_BE = 2.0_double*vdot(nsf*nsf*mat(1,Hrange)%length,data_K,1,data_H,1)
  ! get new H matrix 
  total_energy_1 = total_energy
  call get_H_matrix( .false., fixed_potential, electrons, rhoin, size)
  if(flag_perform_cDFT) then
     call cdft_min (reset_L, fixed_potential, vary_mu, &
          n_CG_L_iterations, tolerance, mu, total_energy)
  else
     ! ** Useful, but not rigorous **
     !new_BE = 2.0_double*vdot(nsf*nsf*mat(1,Hrange)%length,data_K,1,data_H,1)
     !if(inode==ionode) write(io_lun,*) 'Old and new BE are: ',start_BE, new_BE
     !if(reset_L) then
     !   if(abs(new_BE-start_BE)<0.5.AND.abs(new_BE-start_BE)>1.0e-10_double) then
     !      reset_L = .false.
     !      Ltol = 0.05_double * tolerance
     !   else
     !      reset_L = .true.
     !      Ltol = tolerance
     !   end if
     !else
     !   Ltol = 0.05_double * tolerance
     !end if
     Ltol = tolerance
     ! Find minimum density matrix
     call stop_timer(tmr_std_chargescf)
     call FindMinDM(n_CG_L_iterations, vary_mu, Ltol, mu, inode, &
          ionode, reset_L, record)
     call start_timer(tmr_std_chargescf)

     ! If we're using O(N), we only have L, and we need K - if diagonalisation, we have K
     if(.NOT.diagon) call LNV_matrix_multiply(electrons, total_energy_1,  &
          doK, dontM1, dontM2, dontM3, dontM4, dontphi, dontE,0,0,0,0)
     ! Get total energy
     if(flag_check_DFT) then
        call get_energy(total_energy,.false.)
     else
        call get_energy(total_energy)
     endif
  end if

  ! And get that density
  temp_supp_fn = allocate_temp_fn_on_grid(sf)
  call stop_timer(tmr_std_chargescf)      ! This routine is always call within area 5
  call get_electronic_density(rhoout, electrons, supportfns, temp_supp_fn, inode, ionode, size)
  call start_timer(tmr_std_chargescf)  
  call free_temp_fn_on_grid(temp_supp_fn)

  !For DFT energy with charge density constructed by density matrix   --  TM  Nov2007
  if(flag_check_DFT) then
   call get_H_matrix( .false., fixed_potential, electrons, rhoout, size)
   call get_energy(total_energy,flag_check_DFT)
  endif

  return
end subroutine get_new_rho_nospin
!!***


!!****f* EarlySC_module/get_new_rho_spin 
!! PURPOSE 
!!   Gets a new charge density when given an input one for spin
!!   polarised calculations ; we input rho, find the ground-state
!!   density matrix and generate rhoout
!!
!!   Based on EarlySC_module/get_new_rho
!! INPUTS
!! OUTPUT
!! RETURN VALUE
!! AUTHOR
!!   L.Tong
!! CREATION DATE 
!!   2011/09/12
!! MODIFICATION HISTORY
!!   2011/09/13 L.Tong
!!     removed obsolete dependence on number_of_bands
!! SOURCE
!!
subroutine get_new_rho_spin (record, reset_L, fixed_potential, &
     vary_mu, n_CG_L_iterations, tolerance, mu, total_energy, &
     rhoin_up, rhoin_dn, rhoout_up, rhoout_dn, size)
  
  use datatypes
  use logicals
  use mult_module, ONLY: LNV_matrix_multiply
  use DMMin, ONLY: FindMinDM
  use global_module, ONLY: iprint_SC, sf, flag_perform_cDFT
  use H_matrix_module, ONLY: get_H_matrix
  use DiagModule, ONLY: diagon
  use energy, ONLY: get_energy, flag_check_DFT
  use functions_on_grid, ONLY: supportfns, allocate_temp_fn_on_grid, &
       free_temp_fn_on_grid
  use density_module, ONLY: get_electronic_density
  use GenComms, ONLY: inode, ionode
  use cdft_module, ONLY: cdft_min
  
  implicit none

  ! Passed Variables
  logical :: record, reset_L, fixed_potential, vary_mu
  integer :: size, temp_supp_fn
  integer :: n_CG_L_iterations
  real(double) :: tolerance, mu, total_energy
  real(double), dimension(size) :: rhoin_up, rhoin_dn
  real(double), dimension(size) :: rhoout_up, rhoout_dn

  ! Local variables
  real(double) :: electrons_up, electrons_dn, electrons, &
       total_energy_1_up, total_energy_1_dn, start_BE, new_BE, Ltol
  
  call get_H_matrix (.false., fixed_potential, electrons_up, &
       electrons_dn, rhoin_up, rhoin_dn, size)
  
  if (flag_perform_cDFT) then
     call cdft_min (reset_L, fixed_potential, vary_mu, &
          n_CG_L_iterations, tolerance, mu, total_energy)
  else
     Ltol = tolerance
     ! Find minimum density matrix
     call stop_timer (tmr_std_chargescf)
     call FindMinDM (n_CG_L_iterations, vary_mu, Ltol, mu, inode, &
          ionode, reset_L, record)
     call start_timer(tmr_std_chargescf)
     
     ! If we're using O(N), we only have L, and we need K - if
     ! diagonalisation, we have K
     if (.NOT.diagon) then
        call LNV_matrix_multiply (electrons_up, total_energy_1_up, &
             doK, dontM1, dontM2, dontM3, dontM4, dontphi, dontE, 0, &
             0, 0, 0, spin=1)
        call LNV_matrix_multiply (electrons_dn, total_energy_1_dn, &
             doK, dontM1, dontM2, dontM3, dontM4, dontphi, dontE, 0, &
             0, 0, 0, spin=2)
     end if
     ! Get total energy
     if (flag_check_DFT) then
        call get_energy (total_energy, .false.)
     else
        call get_energy (total_energy)
     endif
  end if ! if (flag_perform_cDFT) then

  ! And get the output density
  temp_supp_fn = allocate_temp_fn_on_grid(sf)
  call stop_timer (tmr_std_chargescf) ! This routine is always call within area 5
  call get_electronic_density (rhoout_up, electrons_up, supportfns, &
       temp_supp_fn, inode, ionode, size, spin=1)
  call get_electronic_density (rhoout_dn, electrons_dn, supportfns, &
       temp_supp_fn, inode, ionode, size, spin=2)
  call start_timer(tmr_std_chargescf)  
  call free_temp_fn_on_grid(temp_supp_fn)

  !For DFT energy with charge density constructed by density matrix --
  !TM Nov2007
  if(flag_check_DFT) then
   call get_H_matrix (.false., fixed_potential, electrons_up, &
        electrons_dn, rhoout_up, rhoout_dn, size)
   call get_energy (total_energy, flag_check_DFT)
  endif

  return

end subroutine get_new_rho_spin
!!*****


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
!!  SOURCE
!!
  subroutine mixtwo(len,linear,lambda,rho0, rho1, rhomix)

    use datatypes
    use numbers
    use GenBlas
    use GenComms, ONLY: gsum

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
