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
!!    Improved headers slightly, and converted get_new_rho to check for solution method before
!!    generating K (and not generate it if using diagonalisation)
!!   31/07/2002 dave
!!    Changed matrix_data usage in get_new_rho to use new variable data_M12
!!   2008/02/01 17:47 dave
!!    Changes to write to file not stdout
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

contains

! -----------------------------------------------------------
! Subroutine reduceLambda
! -----------------------------------------------------------

!!****f* EarlySC_module/reduceLambda *
!!
!!  NAME 
!!   reduceLambda
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
!!  SOURCE
!!
  subroutine reduceLambda(R0,R1,lambda_1,thresh,MixLin, reset_L, fixed_potential, vary_mu, n_L_iterations, mx_SCiters, &
       number_of_bands, L_tol, mu, total_energy, rho0,rho1,resid, size)

    use datatypes
    use numbers
    use GenComms, ONLY: cq_abort, inode, ionode

    implicit none

    ! Passed variables
    logical :: MixLin, reset_L, fixed_potential, vary_mu

    integer :: size
    integer :: n_L_iterations, mx_SCiters

    real(double) :: R0, R1, lambda_1, thresh
    real(double) :: number_of_bands, L_tol, mu
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
       R1 = getR2( MixLin, lambda_1, reset_L, fixed_potential, vary_mu, n_L_iterations, number_of_bands, &
            L_tol, mu, total_energy, rho0, rho1,resid, size)
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
!!  SOURCE
!!
  subroutine bracketMin(moved, &
       lambda_a,lambda_b,lambda_c,Rcross_a, Rcross_b, Rcross_c, &
       Ra, Rb, Rc, MixLin, reset_L, fixed_potential, vary_mu, &
       n_L_iterations, mx_SCiters, number_of_bands, L_tol, mu,&
       total_energy, den0,den1, residb, resid0, size)

    use datatypes
    use numbers
    use GenBlas
    use dimens, ONLY: n_my_grid_points
    use GenComms, ONLY: cq_abort, gsum, inode, ionode

    implicit none

    ! Passed variables
    logical :: moved, MixLin ! Flags whether point B has moved
    logical vary_mu, fixed_potential, reset_L

    integer :: size
    integer :: mx_SCiters, n_L_iterations

    real(double) :: lambda_a, lambda_b, lambda_c
    real(double) :: Rcross_a, Rcross_b, Rcross_c, Ra, Rb, Rc
    real(double) :: number_of_bands, L_tol, mu
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
    Rc = getR2( MixLin, lambda_c, reset_L, fixed_potential, vary_mu, n_L_iterations, number_of_bands, L_tol, mu, &
         total_energy, den0, den1,residc, size)
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
             Ru = getR2( MixLin, u, reset_L, fixed_potential, vary_mu, n_L_iterations, number_of_bands, L_tol, mu, &
                  total_energy, den0, den1,residu, size)
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
                Ru = getR2( MixLin, u, reset_L, fixed_potential, vary_mu, n_L_iterations, number_of_bands, L_tol, mu, &
                     total_energy, den0, den1,residu, size)
                if(inode==ionode) write(io_lun,108) n_brak, Ru
                Rcross = dot(n_my_grid_points, residu, 1, resid0, 1)
                call gsum(Rcross)
             endif
          else if((lambda_c - u)*(u-ulim)>zero) then ! Fit between c and ulim
             n_brak = n_brak + 1
             Ru = getR2( MixLin, u, reset_L, fixed_potential, vary_mu, n_L_iterations, number_of_bands, L_tol, mu, &
                  total_energy, den0, den1,residu, size)
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
                Ru = getR2( MixLin, u, reset_L, fixed_potential, vary_mu, n_L_iterations, number_of_bands, L_tol, mu, &
                     total_energy, den0, den1,residu, size)
                Rcross = dot(n_my_grid_points, residu, 1, resid0, 1)
                call gsum(Rcross)
                if(inode==ionode) write(io_lun,108) n_brak, Ru
             endif
          else if((u-ulim)*(ulim-lambda_c)>=zero) then ! Parabolic u to ulim
             u = ulim
             n_brak = n_brak +1
             Ru = getR2( MixLin, u, reset_L, fixed_potential, vary_mu, n_L_iterations, number_of_bands, L_tol, mu, &
                  total_energy, den0, den1,residu, size)
             Rcross = dot(n_my_grid_points, residu, 1, resid0, 1)
             call gsum(Rcross)
             if(inode==ionode) write(io_lun,108) n_brak, Ru
          else ! Reject parabolic fit and do default magnify
             u = lambda_c + GOLD*(lambda_c-lambda_b)
             n_brak = n_brak + 1
             Ru = getR2( MixLin, u, reset_L, fixed_potential, vary_mu, n_L_iterations, number_of_bands, L_tol, mu, &
                  total_energy, den0, den1,residu, size)
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
!!  SOURCE
!!
  subroutine searchMin(lambda_a,lambda_b,lambda_c,Rcross_a, Rcross_c, Ra, Rb, Rc, MixLin, reset_L, &
       fixed_potential, vary_mu, n_L_iterations, mx_SCiters, number_of_bands, L_tol, mu,&
       total_energy, den0,den1, residb, size)

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
    real(double) :: number_of_bands, L_tol, mu
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
       f2 = getR2( MixLin, x2, reset_L, fixed_potential, vary_mu, n_L_iterations, number_of_bands, L_tol, mu, &
            total_energy, den0, den1,resid, size)
       if(inode==ionode) write(io_lun,108) n_s, f2
       Rl = f2
    else
       x2 = lambda_b
       x1 = R*lambda_b + C*lambda_a ! or b-C*(b-a)
       f2 = Rb
       f1 = getR2( MixLin, x1, reset_L, fixed_potential, vary_mu, n_L_iterations, number_of_bands, L_tol, mu, &
            total_energy, den0, den1,resid, size)
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
          f2 = getR2( MixLin, x2, reset_L, fixed_potential, vary_mu, n_L_iterations, number_of_bands, L_tol, mu, &
               total_energy, den0, den1,resid, size)
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
          f1 = getR2( MixLin, x1, reset_L, fixed_potential, vary_mu, n_L_iterations, number_of_bands, L_tol, mu, &
               total_energy, den0, den1,resid, size)
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
!!  SOURCE
!!
  subroutine brentMin(lambda_a,lambda_b,lambda_c,Rcross_a, Rcross_c, Ra, Rb, Rc, MixLin, reset_L, &
       fixed_potential, vary_mu, n_L_iterations, mx_SCiters, number_of_bands, L_tol, mu,&
       total_energy, den0,den1, residb, size)

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
    real(double) :: number_of_bands, L_tol, mu
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
       if(inode==ionode) &
            write(io_lun,*) 'Iteration ',n_s,' in brentMin, x,Rx: ',x,fx
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
       fu = getR2( MixLin, u, reset_L, fixed_potential, vary_mu, n_L_iterations, number_of_bands, L_tol, mu, &
            total_energy, den0, den1,resid, size)
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
  end subroutine brentMin
!!***

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
!!  SOURCE
!!
  real(double) function getR2( linear, lambda, reset_L, fixed_potential, vary_mu, n_L_iterations, &
       number_of_bands, L_tolerance, mu, total_energy, den0, den1,resid, size)

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

    real(double) :: lambda
    real(double) :: number_of_bands, L_tolerance, mu
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
    call get_new_rho( .false., reset_L, fixed_potential, vary_mu, n_L_iterations, number_of_bands, &
         L_tolerance, mu, total_energy, rhoin, rhoout, size)
    resid = zero
    call axpy(size, one, rhoout, 1, resid, 1)
    call axpy(size, -one, rhoin, 1, resid, 1)
    R = dot(size,resid,1,resid,1)
    call gsum(R)
    getR2 = sqrt(grid_point_volume*R)/ne_in_cell
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
!!  SOURCE
!!
subroutine get_new_rho( record, reset_L, fixed_potential, vary_mu, &
     n_CG_L_iterations, number_of_bands, tolerance, mu, total_energy, &
     rhoin, rhoout, size)
  
  use datatypes
  use logicals
  use mult_module, ONLY: LNV_matrix_multiply
  use DMMin, ONLY: FindMinDM
  use global_module, ONLY: iprint_SC, sf
  use H_matrix_module, ONLY: get_H_matrix
  use DiagModule, ONLY: diagon
  use energy, ONLY: get_energy, flag_check_DFT
  use functions_on_grid, ONLY: supportfns, allocate_temp_fn_on_grid, free_temp_fn_on_grid
  use density_module, ONLY: get_electronic_density
  use GenComms, ONLY: inode, ionode
  
  implicit none
  
  ! Passed Variables
  logical :: record, reset_L, fixed_potential, vary_mu
  
  integer :: size, temp_supp_fn
  integer :: n_CG_L_iterations
  
  real(double) :: number_of_bands, tolerance, mu, total_energy
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
  call FindMinDM(n_CG_L_iterations, number_of_bands, vary_mu, Ltol, mu, inode, ionode, reset_L, record)
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
end subroutine get_new_rho
!!***

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
