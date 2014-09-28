! $Id$
! -----------------------------------------------------------
! Module SelfCon_module
! -----------------------------------------------------------
! Code area 5: self-consistency
! -----------------------------------------------------------

!!****h* Conquest/SelfCon_module
!!  NAME
!!   SelfCon_module
!!  PURPOSE
!!   Contains master subroutines used in self-consistency, and various
!!   useful variables.  Both "early" and "late" stage self-consistency
!!   strategies are described in the notes "How to achieve
!!   self-consistency" by MJG to be found in the Conquest
!!   documentation repository, also in Chem. Phys. Lett.  325, 796
!!   (2000) and in pieces in Cquest (my directory NewSC)
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   16/10/00
!!  MODIFICATION HISTORY
!!   16/10/00-06/12/00 by D.R.Bowler
!!    Continuing to implement and improve the routines
!!   18/05/2001 dave
!!    Stripped calls to reduceLambda, getR2 and get_new_rho
!!   08/06/2001 dave
!!    Updated to use GenComms
!!   20/06/2001 dave
!!    Moved fitting of PosTan coefficients to fit_coeff
!!   16:01, 04/02/2003 drb
!!    Small changes: made residuals square roots of dot product (as
!!    they should be)
!!   10:35, 06/03/2003 drb
!!    Added simple linear mixing and one or two other small changes
!!   2006/03/06 05:57 dave
!!    Rewrote calls to entire module to simplify
!!   2006/09/20 17:06 dave
!!    Made parameters user-definable
!!   2008/02/04 08:23 dave
!!    Changed for output to file not stdout
!!   2008/04/02  M. Todorovic
!!    Added atomic charge calculation
!!   2011/09/12 L.Tong
!!    Added A_dn for mixing parameter for spin down component For spin
!!    polarised cases A is used as mixing parameter for spin up
!!    component and A_dn is used as mixing parameter for spin down
!!    component
!!   2012/03/01 L.Tong
!!    Added interface for earlySC and lateSC
!!   2014/09/15 18:30 lat
!!   -fixed call start/stop_timer to timer_module (not timer_stdlocks_module !)
!!  SOURCE
!!
module SelfCon

  use datatypes
  use global_module,          only: io_lun, area_SC
  use timer_module,           only: start_timer, stop_timer, cq_timer
  use timer_stdclocks_module, only: tmr_std_chargescf, tmr_std_allocation

  implicit none

  ! These should all be read or dynamically allocated
  real(double), parameter :: thresh = 2.0_double
  real(double), parameter :: InitialLambda = 1.0_double
  real(double), parameter :: ReduceLimit = 0.5_double
  real(double), parameter :: crit_lin = 0.1_double
  integer :: maxitersSC
  integer :: maxearlySC
  integer :: maxpulaySC
  !integer, parameter :: mx_SCiters = 50
  !integer, parameter :: mx_store_early = 3
  !integer, parameter :: mx_pulay = 5
  integer, allocatable, dimension(:) :: EarlyRecord
  integer, save :: earlyIters
  logical :: earlyL
  logical, parameter :: MixLin = .true.
  integer :: n_exact
  logical :: atomch_output

  logical, save :: flag_Kerker
  logical, save :: flag_wdmetric

  real(double), dimension(2), save :: A ! A(spin) mixing factor for spin
  real(double), save :: q0
  real(double), save :: q1
  logical,      save :: flag_linear_mixing
  real(double), save :: EndLinearMixing

  ! -------------------------------------------------------
  ! RCS ident string for object file id
  ! -------------------------------------------------------
  character(len=80), save, private :: &
       RCSid = "$Id$"
  !!***

contains

  !!****f* SelfCon_module/new_SC_potl *
  !!
  !!  NAME
  !!   new_SC_potl
  !!  USAGE
  !!
  !!  PURPOSE
  !!   This controls the switch from early-stage to late-stage mixing.
  !!   We can also select whether to record the variation of dE vs R for
  !!   the positive tangent code which gives sliding tolerances
  !!
  !!   N.B. We want the charge density in the variable density which is
  !!   stored in density_module.  So on entry, we assume that the incoming
  !!   density is stored in that variable, and the final density will be
  !!   left there
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
  !!    Added ROBODoc header and continued to strip subroutine calls
  !!    Stripped subroutine call
  !!   20/06/2001 dave
  !!    Removed fitting of C and beta to PosTan and called fit_coeff
  !!   10:37, 06/03/2003 drb
  !!    Added call to simple linear mixing
  !!   2007/04/17 09:37 dave
  !!    Changed tolerance for non-self-consistent call to DMM minimisation
  !!   2007/05/01 10:05 dave
  !!    Added converged flag to allow continuation without convergence
  !!   Saturday, 2011/07/30 L.Tong
  !!    Implemented spin polarisation
  !!   2011/09/13 L.Tong
  !!    Removed obsolete dependence on number_of_bands
  !!   2011/09/19 L.Tong
  !!    Only PulayMixSC_spin works correctly for spin polarised
  !!    calculations, and hence edited the subroutine to reflect this
  !!   2012/03/21 L.Tong
  !!   - Major rewrite of spin polarised implementation
  !!   - Removed input parameter real(double) mu
  !!  SOURCE
  !!
  subroutine new_SC_potl(record, self_tol, reset_L, fixed_potential, &
                         vary_mu, n_L_iterations, L_tol, total_energy)

    use datatypes
    use PosTan,          only: PulayC, PulayBeta, SCC, SCBeta, pos_tan, &
                               max_iters, SCE, SCR, fit_coeff
    use numbers
    use global_module,   only: iprint_SC, flag_self_consistent,         &
                               flag_SCconverged, IPRINT_TIME_THRES2,    &
                               nspin, spin_factor
    use H_matrix_module, only: get_H_matrix
    use DMMin,           only: FindMinDM
    use energy,          only: get_energy
    use GenComms,        only: inode, ionode, cq_abort
    use dimens,          only: n_my_grid_points
    use density_module,  only: density
    use maxima_module,   only: maxngrid
    use timer_module

    implicit none

    ! Passed variables
    logical      :: record   ! Flags whether to record dE vs R
    logical      :: vary_mu, fixed_potential, reset_L
    integer      :: n_L_iterations
    real(double) :: self_tol
    real(double) :: L_tol
    real(double) :: total_energy

    ! Local variables
    integer        :: mod_early, ndone, i, nkeep, ndelta, stat, spin
    logical        :: done, problem, early
    real(double)   :: SC_tol, DMM_tol, LastE
    type(cq_timer) :: tmr_l_tmp1
    type(cq_timer) :: tmr_std_loc
    real(double), dimension(nspin) :: electrons

    
!****lat<$
    call start_timer(t=tmr_std_loc,who='new_SC_potl',where=5,level=2,echo=.true.)
!****lat>$

    call start_timer(tmr_std_chargescf)

    ! Build H matrix *with NL and KE*
    call get_H_matrix(.true., fixed_potential, electrons, density, &
                      maxngrid)

    ! The H matrix build is already timed on its own, so I leave out
    ! of the SC preliminaries (open to discussion)
    call start_timer(tmr_l_tmp1,WITH_LEVEL)
!!$    if (.not. flag_self_consistent) then
!!$       call stop_timer(tmr_std_chargescf)
!!$       call FindMinDM(n_L_iterations, vary_mu, L_tol, inode, ionode, &
!!$                      reset_L, .false.)
!!$       call start_timer(tmr_std_chargescf)
!!$       call get_energy(total_energy)
!!$       call stop_print_timer(tmr_l_tmp1, "new_SC_potl (except H build)", &
!!$                             IPRINT_TIME_THRES2)
!!$       call stop_timer(tmr_std_chargescf)
!!$       return
!!$    end if
    if (inode == ionode) &
         write (io_lun, &
                fmt='(8x,"Starting self-consistency.  Tolerance: ",e12.5,/)') &
               self_tol
    if (record) then
       if (inode == ionode .and. iprint_SC > 1) &
            write (io_lun, *) 'Original tol: ', L_tol
       SC_tol  = self_tol
       DMM_tol = L_tol
    else
       SC_tol  = max(self_tol, SCC * (self_tol**SCBeta))
       DMM_tol = max(L_tol, PulayC*((0.1_double * self_tol)**PulayBeta))
    end if
!!$
!!$
!!$
!!$
    ndone   =  0
    done    = .false.
    problem = .false.
    early   = .false.
    flag_SCconverged = .true.
    ! Check on whether we need to do early iterations
    if (.not. allocated(EarlyRecord)) then
       allocate(EarlyRecord(maxearlySC), STAT=stat)
       if (stat /= 0) &
            call cq_abort("Error allocating EarlyRecord in SCpotl: ", maxearlySC)
       EarlyRecord = 0
    end if
    do i = 1, maxearlySC
       if (EarlyRecord(i) == 0) early = .true.
       if (inode == ionode .and. iprint_SC > 1) &
            write(io_lun,*) 'early: ', i, EarlyRecord(i)
    end do
    call stop_print_timer(tmr_l_tmp1, "SCF preliminaries (except H build)", &
                          IPRINT_TIME_THRES2)
!!$
!!$
!!$
!!$
!!$
!!$
!!$
!!$
    ! Loop until self-consistent
    do while (.not. done)
       call start_timer(tmr_l_tmp1, WITH_LEVEL)
       if (flag_linear_mixing) then
          !call LinearMixSC (done,ndone, SC_tol, call
          !call PulayMixSCA (done,ndone, SC_tol, reset_L, fixed_potential,
          !  vary_mu, n_L_iterations, DMM_tol, mu, total_energy,
          !   density, maxngrid)

          call PulayMixSC_spin(done, ndone, SC_tol, reset_L,          &
                               fixed_potential, vary_mu,              &
                               n_L_iterations, DMM_tol, total_energy, &
                               density, maxngrid)

       else if (early .or. problem) then ! Early stage strategy
          earlyL  = .false.
          reset_L = .true.
          if (inode == ionode .and. iprint_SC > 0) &
               write(io_lun,*) '********** EarlySC **********'
          if (inode == ionode .and. problem) &
               write(io_lun,*) 'Problem restart'

          call earlySC(record,done,earlyL,ndone, SC_tol, reset_L, &
                       fixed_potential, vary_mu, n_L_iterations, &
                       DMM_tol, total_energy, density, maxngrid)

          ! Check for early/late stage switching and update variables
          mod_early = 1 + mod(earlyIters,maxearlySC)
          if (earlyL) EarlyRecord(mod_early) = 1
       end if ! (flag_linear_mixing)

       call stop_print_timer(tmr_l_tmp1, "SCF iteration - Early", &
                             IPRINT_TIME_THRES2)

       if (.not. done) then ! Late stage strategy
          call start_timer(tmr_l_tmp1, WITH_LEVEL)
          reset_L = .true. !.false.
          if (inode == ionode .and. iprint_SC > 0) &
               write (io_lun, *) '********** LateSC **********'

          call lateSC(record, done, ndone, SC_tol, reset_L,     &
                      fixed_potential, vary_mu, n_L_iterations, &
                      DMM_tol, total_energy, density, maxngrid)

          problem = .true.  ! This catches returns from lateSC when
          !  not done
          call stop_print_timer(tmr_l_tmp1, "SCF iteration - Late", &
                                IPRINT_TIME_THRES2)
       end if

       ! Now decide whether we need to bother going to early
       early = .false.
       do i = 1, maxearlySC
          if (EarlyRecord(i) == 0) early = .true.
          if (inode == ionode .and. iprint_SC > 1) &
               write (io_lun, *) 'early: ', i, EarlyRecord(i)
       end do
    end do ! while
!!$
!!$
!!$
!!$
!!$
!!$
!!$
!!$
    call start_timer(tmr_l_tmp1, WITH_LEVEL)

    if (record) then ! Fit the C and beta coefficients
       if (inode == ionode .and. iprint_SC > 1) then
          write (io_lun, *) '  List of residuals and energies'
          if (ndone > 0) then
             do i = 1, ndone
                write (io_lun, 7) i, SCR(i), SCE(i)
             end do
          end if
       end if
       call fit_coeff(SCC, SCBeta, SCE, SCR, ndone)
       if (inode == ionode .and. iprint_SC > 1) &
            write (io_lun, 6) SCC, SCBeta
    end if

    call stop_print_timer(tmr_l_tmp1, "finishing SCF (fitting coefficients)", &
                          IPRINT_TIME_THRES2)
    call stop_timer(tmr_std_chargescf)

!****lat<$
    call stop_timer(t=tmr_std_loc,who='new_SC_potl',echo=.true.)
!****lat>$

    return

6   format(8x,'dE to dR parameters - C: ',f15.8,' beta: ',f15.8)
7   format(8x,i4,2f15.8)

  end subroutine new_SC_potl
  !!***


  !!****f* SelfCon_module/earlySC *
  !!
  !!  NAME
  !!   earlySC
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Does the early-stage charge mixing.  This takes a pessimistic view of
  !!   the procedure, and simply performs successive line minimisations until
  !!   the change in the residual is linear in the change in charge (to a
  !!   specified tolerance) before switching to late-stage mixing.
  !!  INPUTS
  !!
  !!
  !!  USES
  !!   EarlySC_module
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!
  !!  MODIFICATION HISTORY
  !!   18/05/2001 dave
  !!    Reduced subroutine calls to get_new_rho, reduceLambda, bracketMin
  !!   08/06/2001 dave
  !!    Changed to use GenComms and gsum
  !!   10:36, 06/03/2003 drb
  !!    Commented out sqrt for R0 as it seems incompatible with EarlySC definitions
  !!   15:01, 02/05/2005 dave
  !!    Changed definition of residual in line with SelfConsistency notes
  !!   2011/09/13 L.Tong
  !!    Removed absolete dependence on number_of_bands
  !!   2012/03/01 L.Tong
  !!    Renamed to earlySC_nospin
  !!   2012/03/21 L.Tong
  !!   - Merged subroutines earlySC_nospin and earlySC_spin into one subroutine
  !!   - Major rewrite of spin implementation
  !!   - Removed redundant input parameter mu
  !!   - Corrected integration measure and normalisation
  !!  SOURCE
  !!
  subroutine earlySC(record, done, EarlyLin, ndone, self_tol, reset_L,&
                     fixed_potential, vary_mu, n_L_iterations, L_tol, &
                     total_energy, rho, size)

    use datatypes
    use numbers
    use GenBlas
    use PosTan
    use EarlySCMod
    use GenComms,      only: cq_abort, gsum, my_barrier, inode, ionode
    use dimens,        only: n_my_grid_points, grid_point_volume
    use global_module, only: ne_in_cell, area_SC, iprint_SC, &
                             flag_continue_on_SC_fail,       &
                             flag_SCconverged, nspin, spin_factor
    use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use maxima_module, only: maxngrid

    implicit none

    ! Passed variables
    logical      :: record, done, EarlyLin
    logical      :: vary_mu, fixed_potential, reset_L
    integer      :: ndone, size
    integer      :: n_L_iterations
    real(double) :: self_tol
    real(double) :: L_tol
    real(double) :: total_energy
    real(double), dimension(:,:) :: rho

    ! Local variables
    integer :: n_iters, irc, stat, spin
    logical :: linear, left, moved, init_reset, Lrec
    real(double) :: R0, R1, Rcross, Rcross_a, Rcross_b, Rcross_c
    real(double) :: Ra, Rb, Rc, E0, E1, dE
    real(double) :: lambda_1, lambda_a, lambda_b, lambda_c
    real(double) :: zeta_exact, ratio, qatio, lambda_up
    real(double) :: pred_Rb, Rup, Rcrossup, zeta_num, zeta
    real(double), dimension(:,:), allocatable :: resid0, rho1, residb

    type(cq_timer)    :: tmr_std_loc
    
!****lat<$
    call start_timer(t=tmr_std_loc,who='earlySC',&
         where=5,level=3,echo=.true.)
!****lat>$

    allocate(resid0(maxngrid,nspin), rho1(maxngrid,nspin), &
             residb(maxngrid,nspin), STAT=stat)
    if (stat /= 0) call cq_abort("earlySC: Error alloc mem: ", maxngrid, nspin)
    call reg_alloc_mem(area_SC, 3*maxngrid*nspin, type_dbl)

    init_reset = reset_L
    call my_barrier()
    ! Test to check for available iterations
    n_iters = ndone
    if (n_iters >= maxitersSC) then
       if (.not. flag_continue_on_SC_fail) &
            call cq_abort('earlySC: Too many self-consisteny iterations: ', &
                          n_iters, maxitersSC)
       flag_SCconverged = .false.
       done = .true.
       return
    end if
    ! Decide on whether or not to record dE vs R for L min
    if (record) then
       Lrec = .true.
    else
       Lrec = .false.
    endif
    ! Compute residual of initial density
    call get_new_rho(Lrec, reset_L, fixed_potential, vary_mu,        &
                     n_L_iterations, L_tol, total_energy, rho, rho1, &
                     maxngrid)
    Lrec = .false.
    R0 = zero
    do spin = 1, nspin
       resid0(1:n_my_grid_points,spin) = &
            rho1(1:n_my_grid_points,spin) - rho(1:n_my_grid_points,spin)
       R0 = R0 + spin_factor * &
            dot(n_my_grid_points, resid0(:,spin), 1, resid0(:,spin), 1)
    end do
    ! cross term
    R0 = R0 + two * &
         dot(n_my_grid_points, resid0(:,1), 1, resid0(:,nspin), 1)
    call gsum(R0)
    R0 = sqrt(grid_point_volume * R0) / ne_in_cell

    if (inode == ionode) write (io_lun, *) 'In EarlySC, R0 is ', R0
    if (R0 < self_tol) then
       done = .true.
       if (inode == ionode) write (io_lun, *) 'Done ! Self-consistent'
    end if
    linear = .false.
    lambda_1 = InitialLambda

    ! Now loop until we achieve a reasonable reduction or linearity
    ! -------------------------------------------------------------
    ! START OF MAIN LOOP
    ! -------------------------------------------------------------

    E0 = total_energy
    dE = L_tol

    do while ((.not. done) .and. (.not. linear) .and. (n_iters < maxitersSC))
       if (init_reset) reset_L = .true.
       n_iters = n_iters + 1
       ! Go to initial lambda and find R

       if (inode == ionode) &
            write (io_lun, *) '********** Early iter ', n_iters

       R1 = getR2(MixLin, lambda_1, reset_L, fixed_potential, vary_mu,&
                  n_L_iterations, L_tol, total_energy, rho, rho1,     &
                  residb, maxngrid)

       Rcross = zero
       do spin = 1, nspin
          Rcross = Rcross + spin_factor * &
                   dot(n_my_grid_points, residb(:,spin), 1, resid0(:,spin), 1)
       end do
       Rcross = Rcross + dot(n_my_grid_points, residb(:,1), 1, resid0(:,nspin), 1)
       Rcross = Rcross + dot(n_my_grid_points, residb(:,nspin), 1, resid0(:,1), 1)
       call gsum(Rcross)
       Rcross = sqrt(grid_point_volume * Rcross) / ne_in_cell

       if (R0 < self_tol) then
          done = .true.
          if (inode == ionode) write (io_lun, *) 'Done ! Self-consistent'
          exit
       end if

       ! Test for acceptable reduction
       if (inode == ionode) write (io_lun, *) 'In EarlySC, R1 is ', R1
       if (R1 <= thresh * R0) then
          if (inode == ionode) write (io_lun, *) 'No need to reduce'
       else
          if (inode == ionode) write (io_lun, *) 'Reducing lambda'
          call reduceLambda(R0, R1, lambda_1, thresh, MixLin, reset_L,&
                            fixed_potential, vary_mu, n_L_iterations, &
                            maxitersSC, L_tol, total_energy, rho,     &
                            rho1, residb, maxngrid)
          Rcross = zero
          do spin = 1, nspin
             Rcross = Rcross + spin_factor * &
                      dot(n_my_grid_points, residb(:,spin), 1, resid0(:,spin), 1)
          end do
          Rcross = Rcross + &
                   dot(n_my_grid_points, residb(:,1), 1, resid0(:,nspin), 1)
          Rcross = Rcross + &
                   dot(n_my_grid_points, residb(:,nspin), 1, resid0(:,1), 1)
          call gsum(Rcross)
          Rcross = sqrt(grid_point_volume * Rcross) / ne_in_cell
       end if ! Preliminary reduction of lambda

       ! At this point, we have two points: 0 and R0, lambda_1 and R1
       ! Now search for a good bracket
       if (R1 > R0) then  ! Take (a,b) = (lam,0)
          left = .false.
          lambda_a = lambda_1
          Ra = R1
          Rcross_a = Rcross
          lambda_b = zero
          Rb = R0
          Rcross_b = R0
          do spin = 1, nspin
             call copy(n_my_grid_points, resid0(:,spin), 1, residb(:,spin), 1)
          end do
       else            ! Take (a,b) = (0,lam)
          left = .true.
          lambda_a = zero
          Ra = R0
          Rcross_a = R0
          lambda_b = lambda_1
          Rb = R1
          Rcross_b = Rcross
          moved = .true.
       end if

       call bracketMin(moved, lambda_a, lambda_b, lambda_c, Rcross_a, &
                       Rcross_b, Rcross_c, Ra, Rb, Rc, MixLin,        &
                       reset_L, fixed_potential, vary_mu,             &
                       n_L_iterations, maxitersSC, L_tol,             &
                       total_energy, rho, rho1, residb, resid0,       &
                       maxngrid)

       ! Test for acceptability - otherwise do line search
       reset_L = .false.
       if (inode == ionode) &
            write (io_lun, *) 'a, b, c: ', &
                  lambda_a, lambda_b, lambda_c, Ra, Rb, Rc
       if (moved .and. Rb < ReduceLimit * R0) then ! Accept
          if (inode == ionode) &
               write (io_lun, *) 'Line search accepted after bracketing', &
                     R0, Rb
       else ! Line search
          if (.not. moved) then
             if (inode == ionode) &
                  write (io_lun, *) 'Search because intern point not moved'
          else
             if (inode == ionode) &
                  write (io_lun, *) 'Search because R reduction too small'
          end if
          ! Do golden section or brent search
          ! call searchMin(lambda_a, lambda_b, lambda_c, Rcross_a, Rcross_c, &
          call brentMin(lambda_a, lambda_b, lambda_c, Rcross_a,     &
                        Rcross_c, Ra, Rb, Rc, MixLin, reset_L,      &
                        fixed_potential, vary_mu, n_L_iterations,   &
                        maxitersSC, L_tol, total_energy, rho, rho1, &
                        residb, maxngrid)
       end if ! end if(reduced)

       if (abs(lambda_a) > abs(lambda_c)) then
          Rcrossup = Rcross_a
          lambda_up = lambda_a
          Rup = Ra
       else
          Rcrossup = Rcross_c
          lambda_up = lambda_c
          Rup = Rc
       end if

       E1 = total_energy
       dE = E1 - E0
       E0 = E1

       ! Now test for linearity
       zeta_exact = abs (Rb - R0)
       ratio = lambda_b / lambda_up
       qatio = one - ratio
       pred_Rb = qatio * qatio * R0 + ratio * ratio * Rup + two * &
                 ratio * qatio * Rcrossup
       zeta_num = abs(Rb - pred_Rb)
       zeta = zeta_num / zeta_exact
       if (inode == ionode) write (io_lun, 131) zeta
       if (zeta < crit_lin) then
          linear = .true.
          if (inode == ionode) write (io_lun, *) 'Linearity fulfilled'
       end if

       ! And record residual and total energy if necessary
       if (record) then
          SCE(n_iters) = total_energy
          SCR(n_iters) = Rb
       end if

       do spin = 1, nspin
          resid0(1:n_my_grid_points,spin) = residb(1:n_my_grid_points,spin)
          call mixtwo(n_my_grid_points, MixLin, lambda_b, rho(:,spin),&
                      rho1(:,spin), residb(:,spin))
          ! rho(:,spin) = zero
          rho(1:n_my_grid_points,spin) = residb(1:n_my_grid_points,spin)
          ! rho1(:,spin) = zero
          rho1(1:n_my_grid_points,spin) = &
               rho(1:n_my_grid_points,spin) + resid0(1:n_my_grid_points,spin)
       end do

       lambda_1 = lambda_b

       R0 = zero
       do spin = 1, nspin
          R0 = R0 + spin_factor * &
               dot(n_my_grid_points, resid0(:,spin), 1, resid0(:,spin), 1)
       end do
       ! cross term
       R0 = R0 + two * &
            dot(n_my_grid_points, resid0(:,1), 1, resid0(:,nspin), 1)
       call gsum(R0)
       R0 = sqrt(grid_point_volume * R0) / ne_in_cell

       if (inode == ionode) write (io_lun, *) 'In EarlySC_spin, R0 is ', R0
       if (R0 < self_tol) then
          done = .true.
          if (inode == ionode) write (io_lun, *) 'Done ! Self-consistent'
       end if

    end do ! end do while(.NOT.linear.AND..NOT.done)

    ! -------------------------------------------------------------
    ! END OF MAIN LOOP
    ! -------------------------------------------------------------
    ! This helps us to decide when to switch from early to late automatically
    if (linear .and. n_iters == 1) EarlyLin = .true.

    ndone = n_iters

    deallocate(resid0, rho1, residb, STAT=stat)
    if (stat /= 0) call cq_abort("earlySC: Error dealloc mem")
    call reg_dealloc_mem(area_SC, 3*maxngrid*nspin, type_dbl)

!****lat<$
    !call stop_timer(t=tmr_std_loc,who='earlySC',echo=.true.)
!****lat>$

    return

131 format('Linearity monitor for this line search: ',e15.6)

  end subroutine earlySC
  !!***


  !!****f* SelfCon_module/lateSC *
  !!
  !!  NAME
  !!   lateSC
  !!  USAGE
  !!
  !!  PURPOSE
  !!   This implements late-stage mixing (using the GR-Pulay algorithm as
  !!   described in Chem. Phys. Lett. 325, 796 (2000) and embellished according
  !!   to the Conquest notes mentioned above).
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
  !!    Stripped subroutine calls to get_new_rho
  !!   08/06/2001 dave
  !!    Changed to use GenComms and gsum
  !!   15:02, 02/05/2005 dave
  !!    Changed definition of residual in line with SelfConsistency notes
  !!   10:34, 13/02/2006 drb
  !!    Removed unnecessary reference to data_H, data_K
  !!   2011/09/13 L.Tong
  !!    Removed absolete dependence on number_of_bands
  !!   2012/03/01 L.Tong
  !!    Renamed to lateSC_nospin
  !!   2012/03/21 L.Tong
  !!   - Merged subroutines lateSC_nospin and lateSC_spin and renamed
  !!     back to lateSC
  !!   - Major rewrite of spin implementation
  !!   - Removed redundant input parameter real(double) mu
  !!  SOURCE
  !!
  subroutine lateSC(record, done, ndone, self_tol, reset_L,          &
                    fixed_potential, vary_mu, n_L_iterations, L_tol, &
                    total_energy, rho, size)

    use datatypes
    use numbers
    use PosTan
    use Pulay,          only: DoPulay
    use GenBlas,        only: dot, rsum
    use dimens,         only: n_my_grid_points, grid_point_volume
    use EarlySCMod,     only: get_new_rho
    use GenComms,       only: gsum, cq_abort, inode, ionode
    use hartree_module, only: kerker
    use global_module,  only: ne_in_cell, area_SC,      &
                              flag_continue_on_SC_fail, &
                              flag_SCconverged,         &
                              flag_fix_spin_population, &
                              ne_spin_in_cell, nspin, spin_factor
    use maxima_module,  only: maxngrid
    use memory_module,  only: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    ! Passed variables
    logical :: record, done
    logical :: vary_mu, fixed_potential, reset_L

    integer :: size
    integer :: ndone
    integer :: n_L_iterations

    real(double) :: self_tol
    real(double) :: L_tol
    real(double) :: total_energy
    real(double), dimension(:,:) :: rho

    ! Local variables
    logical      :: linear
    integer      :: n_iters, n_pulay, npmod, pul_mx, i, j, stat, spin
    real(double) :: R0, R1, E0, E1, dE, tmp_tot
    real(double), dimension(:,:,:), allocatable :: rho_pul, resid_pul
    real(double), dimension(:,:),   allocatable :: rho1
    real(double), dimension(maxpulaySC,maxpulaySC,nspin) :: Aij
    real(double), dimension(maxpulaySC,nspin) :: alph
    real(double), dimension(nspin) :: R, tmp

    type(cq_timer)    :: tmr_std_loc

!****lat<$
    call start_timer(t=tmr_std_loc,who='lateSC',where=5,level=3,echo=.true.)
!****lat>$

    allocate(rho_pul(maxngrid,maxpulaySC,nspin),   &
             resid_pul(maxngrid,maxpulaySC,nspin), &
             rho1(maxngrid,nspin), STAT=stat)
    if (stat /= 0) &
         call cq_abort("late_SC: Error alloc mem: ", maxngrid, maxpulaySC)
    call reg_alloc_mem(area_SC, (2*maxpulaySC+1)*nspin*maxngrid, type_dbl)

    done = .false.
    linear = .true.
    n_iters = ndone
    if (n_iters >= maxitersSC) then
       if (.not. flag_continue_on_SC_fail) &
            call cq_abort('lateSC: too many iterations: ', &
                          n_iters, maxitersSC)
       flag_SCconverged = .false.
       done = .true.
       return
    end if

    ! Compute residual of initial density
    do spin = 1, nspin
       rho_pul(1:n_my_grid_points,1,spin) = rho(1:n_my_grid_points,spin)
    end do

    call get_new_rho(.false., reset_L, fixed_potential, vary_mu,     &
                     n_L_iterations, L_tol, total_energy, rho, rho1, &
                     maxngrid)
    E0 = total_energy
    R0 = zero
    do spin = 1, nspin
       resid_pul(1:n_my_grid_points,1,spin) = &
            rho1(1:n_my_grid_points,spin) - rho(1:n_my_grid_points,spin)
       call kerker(resid_pul(:,1,spin), maxngrid, q0)
       R0 = R0 + spin_factor * &
            dot(n_my_grid_points, resid_pul(:,1,spin), 1, resid_pul(:,1,spin), 1)
    end do
    ! cross terms
    R0 = R0 + two * &
         dot(n_my_grid_points, resid_pul(:,1,1), 1, resid_pul(:,1,nspin), 1)
    call gsum(R0)
    R0 = sqrt(grid_point_volume * R0) / ne_in_cell
    if (inode == ionode) write (io_lun, *) 'Residual is ', R0

    ! Old output becomes new input
    do spin = 1, nspin
       rho(1:n_my_grid_points,spin) = &
            rho(1:n_my_grid_points,spin) + resid_pul(1:n_my_grid_points,1,spin)
    end do

    ! normalise rho
    do spin = 1, nspin
       tmp(spin) = grid_point_volume * rsum(n_my_grid_points, rho(:,spin), 1)
       call gsum(tmp(spin))
    end do
    if (nspin == 1 .or. flag_fix_spin_population) then
       do spin = 1, nspin
          rho(1:n_my_grid_points,spin) = &
               ne_spin_in_cell(spin) * rho(1:n_my_grid_points,spin) / &
               tmp(spin)
       end do
    else
       tmp_tot = spin_factor * sum(tmp)
       do spin = 1, nspin
          rho(1:n_my_grid_points,spin) = &
               ne_in_cell * rho(1:n_my_grid_points,spin) / tmp_tot
       end do
    end if

    ! Now loop until we achieve SC or something goes wrong
    n_pulay = 0
    do while ((.not. done) .and. (linear) .and. (n_iters < maxitersSC))

       if (R0 < one) reset_L = .false.

       n_iters = n_iters + 1

       if (inode == ionode) &
            write (io_lun, *) '********** Late iter ', n_iters

       n_pulay = n_pulay + 1

       ! Storage for pulay charges/residuals
       npmod = mod(n_pulay, maxpulaySC) + 1
       pul_mx = min(n_pulay + 1, maxpulaySC)
       if (inode == ionode) &
            write (io_lun, *) 'npmod, pul_mx: ', npmod, pul_mx

       ! For the present output, find the residual (i.e. R_n^\prime)
       call get_new_rho(.false., reset_L, fixed_potential, vary_mu, &
                        n_L_iterations, L_tol, total_energy, rho,   &
                        rho1, maxngrid)

       do spin = 1, nspin
          rho_pul(1:n_my_grid_points,npmod,spin) = rho(1:n_my_grid_points,spin)
          resid_pul(1:n_my_grid_points,npmod,spin) = &
               rho1(1:n_my_grid_points,spin) - rho(1:n_my_grid_points,spin)
          call kerker(resid_pul(:,npmod,spin), maxngrid, q0)
       end do

       ! Do the pulay-style mixing
       ! Form the A matrix (=<R_i|R_j>)
       do spin = 1, nspin
          do i = 1, pul_mx
             ! diagonal elements of Aij
             R(spin) = dot(n_my_grid_points, resid_pul(:,i,spin), 1, &
                           resid_pul(:,j,spin), 1)
             call gsum(R(spin))
             Aij(i,i,spin) = R(spin)
             ! Aij is symmetric
             if (i > 1) then
                do j = 1, i - 1
                   R(spin) = dot(n_my_grid_points, resid_pul(:,i,spin), &
                                 1, resid_pul(:,j,spin), 1)
                   call gsum(R(spin))
                   Aij(j,i,spin) = R(spin)
                   Aij(i,j,spin) = R(spin)
                end do ! j
             end if ! (i > 1)
          end do ! i
       end do ! spin

       call DoPulay(npmod, Aij, alph, pul_mx, maxpulaySC, inode, &
                    ionode)

       do spin = 1, nspin
          if (inode == ionode) &
               write (io_lun, *) 'alph (spin=', spin, '): ', &
                                 alph(1:pul_mx,spin)
       end do

       ! Build new input density - we could do this wrap around, but
       ! will need rho_up/dn anyway, so might as well use it.
       rho = zero
       do spin = 1, nspin
          do i = 1, pul_mx
             rho(1:n_my_grid_points,spin) = &
                  rho(1:n_my_grid_points,spin) + &
                  alph(i,spin) * rho_pul(1:n_my_grid_points,i,spin)
          end do
       end do

       ! normalise rho
       do spin = 1, nspin
          tmp(spin) = grid_point_volume * rsum(n_my_grid_points, rho(:,spin), 1)
          call gsum(tmp(spin))
       end do
       if (nspin == 1 .or. flag_fix_spin_population) then
          do spin = 1, nspin
             rho(1:n_my_grid_points,spin) = &
                  ne_spin_in_cell(spin) * rho(1:n_my_grid_points,spin) / &
                  tmp(spin)
          end do
       else
          tmp_tot = spin_factor * sum(tmp)
          do spin = 1, nspin
             rho(1:n_my_grid_points,spin) = &
                  ne_in_cell * rho(1:n_my_grid_points,spin) / tmp_tot
          end do
       end if

       do spin = 1, nspin
          rho_pul(1:n_my_grid_points,npmod,spin) = rho(1:n_my_grid_points,spin)
       end do

       ! Generate the residual either exactly or by extrapolation
       if (mod(n_iters, n_exact) == 0) then

          if (inode == ionode) write (io_lun, *) 'Generating rho exactly'
          do spin = 1, nspin
             resid_pul(1:n_my_grid_points,npmod,spin) = &
                  alph(npmod,spin) * resid_pul(1:n_my_grid_points,npmod,spin)
          end do
          if (pul_mx == maxpulaySC) then
             do spin = 1, nspin
                do i = npmod + 1, npmod + maxpulaySC - 1
                   j = mod(i, maxpulaySC)
                   if (j == 0) j = maxpulaySC
                   resid_pul(1:n_my_grid_points,npmod,spin) =      &
                        resid_pul(1:n_my_grid_points,npmod,spin) + &
                        alph(j,spin) * resid_pul(1:n_my_grid_points,j,spin)
                end do
             end do
          else
             do i = 1, npmod - 1
                do spin = 1, nspin
                   resid_pul(1:n_my_grid_points,npmod,spin) = &
                        resid_pul(1:n_my_grid_points,npmod,spin) + &
                        alph(i,spin) * resid_pul(1:n_my_grid_points,i,spin)
                end do
             end do
          end if

          R1 = zero
          do spin = 1, nspin
             R1 = R1 + spin_factor * &
                  dot(n_my_grid_points, resid_pul(:,npmod,spin), 1, &
                      resid_pul(:,npmod,spin), 1)
          end do
          ! cross term
          ! R1 = R1 + two * &
          !      dot(n_my_grid_points, resid_pul(:,npmod,1), 1, &
          !             resid_pul(:,npmod,nspin), 1)
          call gsum(R1)
          R1 = sqrt(grid_point_volume * R1) / ne_in_cell
          if (inode == ionode) &
               write (io_lun, fmt='(8x,"Predicted residual is ",f20.12)') R1

          call get_new_rho(.false., reset_L, fixed_potential, vary_mu,&
                           n_L_iterations, L_tol, total_energy, rho,  &
                           rho1, maxngrid)

          R1 = zero
          do spin = 1, nspin
             resid_pul(1:n_my_grid_points,npmod,spin) = &
                  rho1(1:n_my_grid_points,spin) - rho(1:n_my_grid_points,spin)
             call kerker(resid_pul(:,npmod,spin), maxngrid, q0)
             R1 = R1 + spin_factor * &
                  dot (n_my_grid_points, resid_pul(:,npmod,spin), 1, &
                       resid_pul(:,npmod,spin), 1)
          end do
          ! cross term
          ! R1 = R1 + two * &
          !      dot(n_my_grid_points, resid_pul(:,npmod,1), 1, &
          !          resid_pul(:,npmod,nspin), 1)
          call gsum(R1)
          R1 = sqrt(grid_point_volume * R1) / ne_in_cell
          if (inode == ionode) &
               write (io_lun, fmt='(8x,"Actual residual is ", f20.12)') R1

          do spin = 1, nspin
             rho(1:n_my_grid_points,spin) = &
                  rho(1:n_my_grid_points,spin) + &
                  resid_pul(1:n_my_grid_points,npmod,spin)
          end do

          ! normalise rho
          do spin = 1, nspin
             tmp(spin) = grid_point_volume * rsum(n_my_grid_points, rho(:,spin), 1)
             call gsum(tmp(spin))
          end do
          if (nspin == 1 .or. flag_fix_spin_population) then
             do spin = 1, nspin
                rho(1:n_my_grid_points,spin) = &
                     ne_spin_in_cell(spin) * rho(1:n_my_grid_points,spin) / &
                     tmp(spin)
             end do
          else
             tmp_tot = spin_factor * sum(tmp)
             do spin = 1, nspin
                rho(1:n_my_grid_points,spin) = &
                     ne_in_cell * rho(1:n_my_grid_points,spin) / tmp_tot
             end do
          end if

       else ! if (mod(n_iters, n_exact) == 0) then

          if (inode == ionode) &
               write (io_lun, *) 'Generating rho by interpolation'
          ! Clever, wrap-around way to do it
          do spin = 1, nspin
             resid_pul(1:n_my_grid_points,npmod,spin) = &
                  alph(npmod,spin) * resid_pul(1:n_my_grid_points,npmod,spin)
          end do
          if (pul_mx == maxpulaySC) then
             do i = npmod + 1, npmod + maxpulaySC - 1
                j = mod(i, maxpulaySC)
                if (j == 0) j = maxpulaySC
                do spin = 1, nspin
                   resid_pul(1:n_my_grid_points,npmod,spin) = &
                        resid_pul(1:n_my_grid_points,npmod,spin) + &
                        alph(j,spin) * resid_pul(1:n_my_grid_points,j,spin)
                end do
             end do
          else
             do spin = 1, nspin
                do i = 1, npmod - 1
                   resid_pul(1:n_my_grid_points,npmod,spin) = &
                        resid_pul(1:n_my_grid_points,npmod,spin) + &
                        alph(i,spin) * resid_pul(1:n_my_grid_points,i,spin)
                end do
             end do
          end if

          ! Output for free if we've extrapolated
          do spin = 1, nspin
             rho(1:n_my_grid_points,spin) = &
                  rho_pul(1:n_my_grid_points,npmod,spin) + &
                  resid_pul(1:n_my_grid_points,npmod,spin)
          end do

          ! normalise rho
          do spin = 1, nspin
             tmp(spin) = grid_point_volume * &
                         rsum(n_my_grid_points, rho(:,spin), 1)
             call gsum(tmp(spin))
          end do
          if (nspin == 1 .or. flag_fix_spin_population) then
             do spin = 1, nspin
                rho(1:n_my_grid_points,spin) = &
                     ne_spin_in_cell(spin) * rho(1:n_my_grid_points,spin) / &
                     tmp(spin)
             end do
          else
             tmp_tot = spin_factor * sum(tmp)
             do spin = 1, nspin
                rho(1:n_my_grid_points,spin) = &
                     ne_in_cell * rho(1:n_my_grid_points,spin) / tmp_tot
             end do
          end if

       end if ! if (mod(n_iters, n_exact) == 0)

       ! Test for (a) convergence and (b) linearity
       E1 = total_energy
       dE = E1 - E0
       E0 = E1
       R1 = zero
       do spin = 1, nspin
          R1 = R1 + spin_factor * &
               dot(n_my_grid_points, resid_pul(:,npmod,spin), 1, &
                   resid_pul(:,npmod,spin), 1)
       end do
       ! R1 = R1 + two * &
       !      dot(n_my_grid_points, resid_pul(:,npmod,1), 1, &
       !          resid_pul(:,npmod,nspin), 1)
       call gsum(R1)
       R1 = sqrt(grid_point_volume * R1) / ne_in_cell

       if (R1 > two * R0 .and. mod(n_iters, n_exact) /= 0) then
          linear = .false.
          if (npmod > 1) then
             do spin = 1, nspin
                rho(1:n_my_grid_points,spin) = &
                     rho_pul(1:n_my_grid_points,npmod-1,spin)
             end do
          else
             do spin = 1, nspin
                rho(1:n_my_grid_points,spin) = &
                     rho_pul(1:n_my_grid_points,pul_mx,spin)
             end do
          end if
          if (inode == ionode) write (io_lun, *) 'PANIC ! Residual increase !'
       end if
       !
       R0 = R1
       !
       if (R0 < self_tol) then
          done = .true.
          if (inode == ionode) write (io_lun, *) 'Done ! Self-consistent'
       end if
       if (inode == ionode) write (io_lun, *) 'Residual is ', R0
       if (record) then
          SCE(n_iters) = total_energy
          SCR(n_iters) = R0
       end if
    end do

    ndone = n_iters

    if (inode == ionode) &
         write(io_lun, *) 'Finishing lateSC after ', ndone, &
                          'iterations with residual of ', R0

    deallocate(rho_pul, resid_pul, rho1, STAT=stat)
    if (stat /= 0) &
         call cq_abort("late_SC: Error dealloc mem")
    call reg_dealloc_mem(area_SC, (2*maxpulaySC+1)*nspin*maxngrid, type_dbl)

!****lat<$
    call stop_timer(t=tmr_std_loc,who='lateSC',echo=.true.)
!****lat>$

    return
  end subroutine lateSC
  !!***


  !!****f* SelfCon_module/LinearMixSC *
  !!
  !!  NAME
  !!   LinearMixSC - simple linear mixing
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Performs simple linear mixing of charge density in Conquest
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   16:52, 2003/03/04
  !!  MODIFICATION HISTORY
  !!   10:11, 12/03/2003 drb
  !!    Made terminating point a parameter
  !!   15:05, 02/05/2005 dave
  !!    Changed definition of residual in line with SelfConsistency notes
  !!   2011/09/13 L.Tong
  !!    Removed obsolete dependence on number_of_bands
  !!   2012/03/21 L.Tong
  !!   - Added spin polarisation
  !!   - Added kerker and wave-dependent metric
  !!   - Removed redundant input parameter real(double) mu
  !!   2013/07/10 11:26 dave
  !!   Bug fix for sum over two components of rho even without spin (and moved rho_total alloc/dealloc)
  !!  SOURCE
  !!
  subroutine LinearMixSC(done, ndone, self_tol, reset_L,           &
                         fixed_potential, vary_mu, n_L_iterations, &
                         L_tol, total_energy, rho, size)

    use datatypes
    use numbers
    use PosTan
    use Pulay,          only: DoPulay
    use GenBlas
    use dimens,         only: n_my_grid_points, grid_point_volume
    use EarlySCMod,     only: get_new_rho
    use GenComms,       only: gsum, cq_abort, inode, ionode
    use io_module,      only: dump_charge
    use hartree_module, only: kerker, kerker_and_wdmetric, wdmetric
    use global_module,  only: ne_in_cell, area_SC, iprint_SC, &
                              flag_continue_on_SC_fail,       &
                              flag_SCconverged, nspin, spin_factor
    use maxima_module,  only: maxngrid
    use memory_module,  only: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    ! Passed variables
    logical      :: done
    logical      :: vary_mu, fixed_potential, reset_L
    integer      :: ndone, size
    integer      :: n_L_iterations
    real(double) :: self_tol
    real(double) :: L_tol
    real(double) :: total_energy
    real(double), dimension(:,:) :: rho

    ! Local variables
    integer      :: n_iters, stat, spin
    real(double) :: R0
    real(double), dimension(:),   allocatable :: rho_tot
    real(double), dimension(:,:), allocatable :: rho1, resid
    real(double), dimension(:,:), allocatable :: resid_cov

    type(cq_timer)    :: tmr_std_loc

!****lat<$
    call start_timer(t=tmr_std_loc,who='LinearMixSC',where=1,level=0,echo=.true.)
!****lat>$

    allocate(rho1(maxngrid,nspin), &
             resid(maxngrid,nspin), STAT=stat)
    if (stat /= 0) call cq_abort("LinearMixSC: Error alloc mem: ", maxngrid, nspin)
    call reg_alloc_mem(area_SC, 2*nspin*maxngrid, type_dbl)

    if (flag_wdmetric) then
       allocate(resid_cov(maxngrid,nspin), STAT=stat)
       if (stat /= 0) &
            call cq_abort("LinearMixSC: failed to allocate resid_cov", stat)
       call reg_alloc_mem(area_SC, nspin * maxngrid, type_dbl)
       resid_cov = zero
    end if

    done = .false.
    n_iters = ndone

    if (n_iters >= maxitersSC) then
       if (.not. flag_continue_on_SC_fail) &
            call cq_abort('LinearMixSC: too many iterations: ', &
                          n_iters, maxitersSC)
       flag_SCconverged = .false.
       done = .true.
       return
    end if

    R0 = 100.0_double
    do while (R0 > EndLinearMixing)

       ! Generate new charge and find residual
       call get_new_rho(.false., reset_L, fixed_potential, vary_mu, &
                        n_L_iterations, L_tol, total_energy, rho,   &
                        rho1, maxngrid)
       do spin = 1, nspin
          resid(1:n_my_grid_points,spin) = &
               rho1(1:n_my_grid_points,spin) - rho(1:n_my_grid_points,spin)
       end do
       if (flag_Kerker) then
          if (flag_wdmetric) then
             R0 = zero
             do spin = 1, nspin
                call kerker_and_wdmetric(resid(:,spin), resid_cov(:,spin), &
                                         maxngrid, q0, q1)
                R0 = R0 + spin_factor * &
                     dot(n_my_grid_points, resid_cov(:,spin), 1, resid(:,spin), 1)
             end do
             ! cross term
             R0 = R0 + dot(n_my_grid_points, resid_cov(:,1), 1, resid(:,nspin), 1)
             R0 = R0 + dot(n_my_grid_points, resid_cov(:,nspin), 1, resid(:,1), 1)
          else
             R0 = zero
             do spin = 1, nspin
                call kerker(resid(:,spin), maxngrid, q0)
                R0 = R0 + spin_factor * &
                     dot(n_my_grid_points, resid(:,spin), 1, resid(:,spin), 1)
             end do
             ! cross term
             R0 = R0 + two * &
                  dot(n_my_grid_points, resid(:,1), 1, resid(:,nspin), 1)
          end if
       else
          if (flag_wdmetric) then
             R0 = zero
             do spin = 1, nspin
                call wdmetric(resid(:,spin), resid_cov(:,spin), maxngrid, q1)
                R0 = R0 + spin_factor * &
                     dot(n_my_grid_points, resid_cov(:,spin), 1, &
                         resid(:,spin), 1)
             end do
             ! cross term
             R0 = R0 + dot(n_my_grid_points, resid_cov(:,1), 1, resid(:,nspin), 1)
             R0 = R0 + dot(n_my_grid_points, resid_cov(:,nspin), 1, resid(:,1), 1)
          else
             R0 = zero
             do spin = 1, nspin
                R0 = R0 + spin_factor * &
                     dot(n_my_grid_points, resid(:,spin), 1, resid(:,spin), 1)
             end do
             ! cross term
             R0 = R0 + two * &
                  dot(n_my_grid_points, resid(:,1), 1, resid(:,nspin), 1)
          end if
       end if

       call gsum(R0)
       R0 = sqrt(grid_point_volume * R0) / ne_in_cell

       if (inode == ionode) write (io_lun, *) 'Residual is ', R0

       ! Old output becomes new input
       !rho = A*rho1 + (1.0_double-A)*rho
       do spin = 1, nspin
          rho(1:n_my_grid_points,spin) = &
               rho(1:n_my_grid_points,spin) + &
               A(spin) * resid(1:n_my_grid_points,spin)
       end do

       n_iters = n_iters + 1

       ! print out charge
       if (iprint_SC > 1) then
          if (nspin == 1) then
             allocate(rho_tot(maxngrid), STAT=stat)
             if (stat /= 0) call cq_abort("LinearMixSC: Error alloc mem: ", maxngrid)
             call reg_alloc_mem(area_SC, maxngrid, type_dbl)
             rho_tot = zero
             rho_tot(:) = spin_factor * rho(:,1)
             call dump_charge(rho_tot, n_my_grid_points, inode, spin=0)
             deallocate(rho_tot, STAT=stat)
             if (stat /= 0) call cq_abort("LinearMixSC: Error dealloc mem")
             call reg_dealloc_mem(area_SC, maxngrid, type_dbl)
          else
             call dump_charge(rho(:,1), n_my_grid_points, inode, spin=1)
             call dump_charge(rho(:,2), n_my_grid_points, inode, spin=2)
          end if
       end if

    end do ! while

    ! deallocate arrays
    if (flag_wdmetric) then
       deallocate(resid_cov, STAT=stat)
       if (stat /= 0) &
            call cq_abort("LinearMixSC: failed to deallocate resid_cov", stat)
       call reg_dealloc_mem(area_SC, nspin * maxngrid, type_dbl)
    end if

    deallocate(rho1, resid, STAT=stat)
    if (stat /= 0) call cq_abort("LinearMixSC: Error dealloc mem")
    call reg_dealloc_mem(area_SC, 2*nspin*maxngrid, type_dbl)

    !    ndone = n_iters

!****lat<$
    call stop_timer(t=tmr_std_loc,who='LinearMixSC',echo=.true.)
!****lat>$

    return
  end subroutine LinearMixSC
  !!***


  ! !!****f* SelfCon_module/PulayMixSC
  ! !! PURPOSE
  ! !!   VASP version
  ! !! INPUTS
  ! !! OUTPUT
  ! !! RETURN VALUE
  ! !! AUTHOR
  ! !!   David Bowler
  ! !! CREATION DATE
  ! !!
  ! !! MODIFICATION HISTORY
  ! !!   2011/09/13 L.Tong
  ! !!     Removed obsolete dependence on number_of_bands
  ! !! SOURCE
  ! !!
  ! subroutine PulayMixSC(done,ndone,self_tol, reset_L, fixed_potential,&
  !      vary_mu, n_L_iterations, L_tol, mu, total_energy,rho,size)
  !   use datatypes
  !   use numbers
  !   use PosTan
  !   use Pulay, only: DoPulay
  !   use GenBlas
  !   use dimens, only: n_my_grid_points, grid_point_volume
  !   use EarlySCMod, only: get_new_rho
  !   use GenComms, only: gsum, cq_abort, inode, ionode
  !   use io_module, only: dump_charge
  !   use hartree_module, only: kerker
  !   use global_module, only: ne_in_cell, area_SC,flag_continue_on_SC_fail, flag_SCconverged
  !   use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl
  !   use maxima_module, only: maxngrid

  !   implicit none

  !   ! Passed variables
  !   logical :: done
  !   logical :: vary_mu, fixed_potential, reset_L

  !   integer :: size
  !   integer :: ndone
  !   integer :: n_L_iterations

  !   real(double) :: self_tol
  !   real(double) :: L_tol, mixing, mu
  !   real(double) :: total_energy
  !   real(double), dimension(size) :: rho

  !   ! Local variables
  !   integer :: n_iters, n_pulay, npmod, pul_mx, i, j, ii, m, next_i, stat
  !   real(double) :: R0, R, R1, E0, E1, dE, stLtol, norm
  !   real(double), allocatable, dimension(:,:) :: rho_pul, delta_rho, delta_R
  !   real(double), allocatable, dimension(:,:) :: resid_pul
  !   real(double), allocatable, dimension(:) :: rho1, resid, Kresid, old_resid
  !   real(double), allocatable, dimension(:,:) :: Aij
  !   real(double), allocatable, dimension(:) :: alph

  !   allocate(rho_pul(maxngrid,maxpulaySC), resid_pul(maxngrid,maxpulaySC),rho1(maxngrid), &
  !        delta_rho(maxngrid,maxpulaySC), delta_R(maxngrid,maxpulaySC), &
  !        resid(maxngrid), Kresid(maxngrid),old_resid(maxngrid), &
  !        Aij(maxpulaySC, maxpulaySC), alph(maxpulaySC), STAT=stat)
  !   if(stat/=0) call cq_abort("Allocation error in PulayMixSC: ",maxngrid, maxpulaySC)
  !   call reg_alloc_mem(area_SC, 4*maxngrid*(1+maxpulaySC)+maxpulaySC*(maxpulaySC+1), type_dbl)
  !   delta_R = zero
  !   delta_rho = zero
  !   done = .false.
  !   n_iters = ndone
  !   if(n_iters>=maxitersSC) then
  !      if(.NOT.flag_continue_on_SC_fail) call cq_abort('lateSC: too many iterations: ',n_iters, maxitersSC)
  !      flag_SCconverged = .false.
  !      done = .true.
  !      return
  !   endif
  !   ! m=1
  !   call get_new_rho(.false., reset_L, fixed_potential, vary_mu, &
  !        n_L_iterations, L_tol, mu, total_energy, rho, rho1, size)
  !   ! Evaluate residual
  !   resid = rho1 - rho
  !   !Kresid = resid
  !   !call kerker(Kresid,maxngrid,q0)
  !   call kerker(resid,maxngrid,q0)
  !   R0 = dot(n_my_grid_points,resid,1,resid,1)
  !   call gsum(R0)
  !   R0 = sqrt(grid_point_volume*R0)/ne_in_cell
  !   if(inode==ionode) write(io_lun,*) 'Residual is ',R0
  !   ! Create new input charge
  !   !rho = rho + A*Kresid
  !   rho = rho + A*resid
  !   ! Store deltarho
  !   delta_rho = zero
  !   call axpy(n_my_grid_points,A,resid,1,delta_rho(:,1),1)
  !   !call axpy(n_my_grid_points,A,Kresid,1,delta_rho(:,1),1)
  !   n_iters = n_iters+1
  !   call dump_charge(rho,size,inode)
  !   do m=2,maxitersSC
  !      ! calculate i (cyclical index for storing history) and pul_mx
  !      i = mod(m-2, maxpulaySC)+1
  !      next_i = mod(i,maxpulaySC)+1
  !      pul_mx = min(m-1, maxpulaySC)
  !      if(inode==ionode) write(io_lun,fmt='(8x,"Pulay iteration ",i4," counters are ",2i4)') m, i, pul_mx
  !      ! Generate new charge and find residual
  !      call get_new_rho(.false., reset_L, fixed_potential, vary_mu, &
  !           n_L_iterations, L_tol, mu, total_energy, rho, rho1, size)
  !      call dump_charge(rho,size,inode)
  !      old_resid = resid
  !      resid = rho1 - rho ! F(m)
  !      call kerker(resid,n_my_grid_points,q0)
  !      R0 = dot(n_my_grid_points,resid,1,resid,1)
  !      call gsum(R0)
  !      R0 = sqrt(grid_point_volume*R0)/ne_in_cell
  !      if(inode==ionode) write(io_lun,*) 'Residual is ',R0
  !      if(R0<self_tol) then
  !         if(inode==ionode) write(io_lun,'(8x,"Reached tolerance")')
  !         done = .true.
  !         call dealloc_PulayMixSC
  !         return
  !      end if
  !      if(R0<EndLinearMixing) then
  !         if(inode==ionode) write(io_lun,'(8x,"Reached transition to LateSC")')
  !         call dealloc_PulayMixSC
  !         return
  !      end if
  !      ! Build deltaR(i)
  !      delta_R(:,i) = zero
  !      call axpy(n_my_grid_points,one,resid,1,delta_R(:,i),1)
  !      call axpy(n_my_grid_points,-one,old_resid,1,delta_R(:,i),1)
  !      !call kerker(resid,n_my_grid_points,q0)
  !      ! Normalise
  !      !norm = dot(n_my_grid_points,delta_R(:,i),1,delta_R(:,i),1)
  !      !call gsum(norm)
  !      !norm = sqrt(norm)
  !      !delta_R(:,i) = delta_R(:,i)/norm
  !      !delta_rho(:,i) = delta_rho(:,i)/norm
  !      !if(inode==ionode) write(io_lun,fmt='(8x,"Norm of deltaF is ",f20.12)') norm
  !      ! now build new rho
  !      do ii=1,pul_mx
  !         do j=1,pul_mx
  !            R = dot(n_my_grid_points, delta_R(:,ii), 1, delta_R(:,j),1)
  !            call gsum(R)
  !            Aij(j,ii) = R
  !            if(ii==j) Aij(j,ii) = Aij(j,ii) !+ 0.0001_double ! w_0 = 0.01
  !         enddo
  !      enddo
  !      if(inode==ionode) write(io_lun,*) 'A is ',Aij
  !      ! Solve to get alphas
  !      call DoPulay(Aij,alph,pul_mx,maxpulaySC,inode,ionode)
  !      if(inode==ionode) write(io_lun,*) 'A is ',Aij
  !      alph = zero
  !      do j=1,pul_mx
  !         R = dot(n_my_grid_points, delta_R(:,j),1,resid(:),1)
  !         call gsum(R)
  !         do ii=1,pul_mx
  !            alph(ii) = alph(ii) + Aij(j,ii)*R
  !         enddo
  !      enddo
  !      if(inode==ionode)write(io_lun,*) 'alph: ',alph(1:pul_mx)
  !      ! Create delta rho - add on new rho
  !      !Kresid = resid
  !      !call kerker(Kresid,maxngrid,q0)
  !      !rho = rho + A*Kresid
  !      !delta_rho(:,next_i) = A*Kresid
  !      rho = rho + A*resid
  !      delta_rho(:,next_i) = A*resid
  !      do ii=1,pul_mx
  !         !Kresid = delta_R(:,ii)
  !         !call kerker(Kresid,maxngrid,q0)
  !         !rho = rho - alph(ii)*(delta_rho(:,ii)+A*Kresid)
  !         !delta_rho(:,next_i) = delta_rho(:,next_i) - alph(ii)*(delta_rho(:,ii)+A*Kresid)
  !         rho = rho - alph(ii)*(delta_rho(:,ii)+A*delta_R(:,ii))
  !         delta_rho(:,next_i) = delta_rho(:,next_i) - alph(ii)*(delta_rho(:,ii)+A*delta_R(:,ii))
  !      end do
  !      n_iters = n_iters+1
  !      call dump_charge(rho,size,inode)
  !   end do
  !   call dealloc_PulayMixSC
  !   !    ndone = n_iters
  !   return

  ! contains

  !   subroutine dealloc_PulayMixSC

  !     deallocate(rho_pul, resid_pul,rho1, delta_rho, delta_R, resid, Kresid,old_resid, Aij, alph, STAT=stat)
  !     if(stat/=0) call cq_abort("Deallocation error in PulayMixSC: ",maxngrid, maxpulaySC)
  !     call reg_dealloc_mem(area_SC, 4*maxngrid*(1+maxpulaySC)+maxpulaySC*(maxpulaySC+1), type_dbl)
  !     return
  !   end subroutine dealloc_PulayMixSC

  ! end subroutine PulayMixSC
  ! !!*****


  ! !!****f* SelfCon_module/PulayMixSCA
  ! !! PURPOSE
  ! !! INPUTS
  ! !! OUTPUT
  ! !! RETURN VALUE
  ! !! AUTHOR
  ! !!   David Bowler
  ! !! CREATION DATE
  ! !!
  ! !! MODIFICATION HISTORY
  ! !!   2011/09/13 L.Tong
  ! !!     Removed obsolete dependence on number_of_bands
  ! !! SOURCE
  ! !!
  ! subroutine PulayMixSCA(done,ndone,self_tol, reset_L, &
  !      fixed_potential, vary_mu, n_L_iterations, L_tol, mu, &
  !      total_energy,rho, size)
  !   use datatypes
  !   use numbers
  !   use PosTan
  !   use Pulay, only: DoPulay
  !   use GenBlas
  !   use dimens, only: n_my_grid_points, grid_point_volume
  !   use EarlySCMod, only: get_new_rho, mixtwo
  !   use GenComms, only: gsum, cq_abort, inode, ionode
  !   use io_module, only: dump_charge, dump_charge2
  !   use hartree_module, only: kerker
  !   use global_module, only: ne_in_cell, iprint_SC, area_SC, flag_continue_on_SC_fail, flag_SCconverged
  !   use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl
  !   use maxima_module, only: maxngrid

  !   implicit none

  !   ! Passed variables
  !   logical :: done
  !   logical :: vary_mu, fixed_potential, reset_L

  !   integer :: size
  !   integer :: ndone
  !   integer :: n_L_iterations

  !   real(double) :: self_tol
  !   real(double) :: L_tol, mixing, mu
  !   real(double) :: total_energy
  !   real(double), dimension(size) :: rho

  !   ! Local variables
  !   integer :: n_iters, n_pulay, npmod, pul_mx, i, j, ii, m, stat, thispulaySC
  !   real(double) :: R0, R1, E0, E1, dE, stLtol, tmp, max_neg, minA, maxA

  !   real(double), allocatable, dimension(:,:) :: rho_pul, R
  !   real(double), allocatable, dimension(:) :: rho1, resid, Kresid
  !   real(double), allocatable, dimension(:,:) :: Aij
  !   real(double), allocatable, dimension(:) :: alph

  !   character(len=11) :: digitstr = "01234567890"

  !   !Reset Pulay Iterations  -- introduced by TM, Nov2007
  !   integer, parameter :: mx_fail = 3
  !   integer :: IterPulayReset=1, icounter_fail=0
  !   logical :: reset_Pulay=.false.
  !   real(double) :: R0_old

  !   allocate(rho_pul(maxngrid,maxpulaySC), R(maxngrid,maxpulaySC),rho1(maxngrid), &
  !        resid(maxngrid), Kresid(maxngrid), Aij(maxpulaySC, maxpulaySC), alph(maxpulaySC), STAT=stat)
  !   if(stat/=0) call cq_abort("Allocation error in PulayMixSC: ",maxngrid, maxpulaySC)
  !   call reg_alloc_mem(area_SC, maxngrid+2*maxngrid*(1+maxpulaySC)+maxpulaySC*(maxpulaySC+1), type_dbl)
  !   rho_pul = zero
  !   R = zero
  !   rho1 = zero
  !   resid = zero
  !   Kresid = zero
  !   Aij = zero
  !   alph = zero
  !   if(inode==ionode) write(io_lun,fmt='(8x,"Starting Pulay mixing, A = ",f6.3," q0= ",f7.4)') A, q0
  !   done = .false.
  !   n_iters = ndone
  !   m=1
  !   if(n_iters>=maxitersSC) then
  !      if(.NOT.flag_continue_on_SC_fail) call cq_abort('lateSC: too many iterations: ',n_iters, maxitersSC)
  !      flag_SCconverged = .false.
  !      done = .true.
  !      return
  !   endif
  !   ! m=1
  !   max_neg = 0.0_double
  !   do j=1,n_my_grid_points
  !      rho_pul(j,1) = rho(j)
  !      if(rho(j)<zero) max_neg = max(max_neg,abs(rho(j)))
  !   end do
  !   if(inode==ionode.AND.max_neg>0.0_double.AND.iprint_SC>3) &
  !        write(io_lun,fmt='(8x,"Init Max negative dens on node ",i5," &
  !        &is ",f8.3)') inode,max_neg
  !   call get_new_rho(.false., reset_L, fixed_potential, vary_mu, &
  !        n_L_iterations, L_tol, mu, total_energy, rho, rho1, size)
  !   ! Evaluate residual
  !   do j=1,n_my_grid_points
  !      resid(j) = rho1(j) - rho(j)
  !   end do
  !   !tmp = grid_point_volume*asum(n_my_grid_points,resid,1)
  !   !write(io_lun,*) 'Sum of resid: ',tmp
  !   R0 = dot(n_my_grid_points,resid,1,resid,1)
  !   call gsum(R0)
  !   R0 = sqrt(grid_point_volume*R0)/ne_in_cell
  !   if(R0<self_tol) then
  !      if(inode==ionode) write(io_lun,'(8x,"Reached self-consistency tolerance")')
  !      done = .true.
  !      call dealloc_PulayMiXSCA
  !      return
  !   end if
  !   if(R0<EndLinearMixing) then
  !      if(inode==ionode) write(io_lun,'(8x,"Reached transition to LateSC")')
  !      call dealloc_PulayMiXSCA
  !      return
  !   end if
  !   !if(inode==ionode.AND.iprint_SC>=0) write(io_lun,*) 'Initial residual is ',R0
  !   if(inode==ionode) write(io_lun,fmt='(8x,"Pulay iteration ",i5," Residual is ",e12.5,/)') m,R0
  !   do j=1,n_my_grid_points
  !      R(j,1) = resid(j)
  !   end do
  !   ! Create new input charge
  !   Kresid = 0.0_double
  !   do j=1,n_my_grid_points
  !      Kresid(j) = resid(j)
  !   end do
  !   call kerker(Kresid,maxngrid,q0)
  !   !tmp = grid_point_volume*asum(n_my_grid_points,Kresid,1)
  !   !write(io_lun,*) 'Sum of Kresid: ',tmp
  !   max_neg = 0.0_double
  !   rho1 = zero
  !   rho1 = rho + Kresid
  !   !tmp = grid_point_volume*asum(n_my_grid_points,rho1,1)
  !   !call gsum(tmp)
  !   !write(io_lun,*) 'Sum of rho1: ',tmp
  !   !rho1 = ne_in_cell*rho1/tmp
  !   !write(io_lun,*) 'Calling mixtwo'
  !   call mixtwo(n_my_grid_points, .true., A, rho, rho1, resid)
  !   !tmp = grid_point_volume*asum(n_my_grid_points,resid,1)
  !   !write(io_lun,*) 'Sum of rho1: ',tmp
  !   rho = resid!*ne_in_cell/tmp
  !   !rho = rho1
  !   do j=1,n_my_grid_points
  !      if(rho(j)<0.0_double) then
  !         if(-rho(j)>max_neg) max_neg = -rho(j)
  !         !   rho(j) = 0.0_double
  !      end if
  !   end do
  !   if(max_neg>0.0_double.AND.iprint_SC>1.AND.inode==ionode) write(io_lun,*) 'First Max negative dens on node ',inode,max_neg
  !   ! Store deltarho
  !   n_iters = n_iters+1
  !   !Reset Pulay Iterations  -- introduced by TM, Nov2007
  !   R0_old = R0
  !   IterPulayReset=1
  !   icounter_fail = 0
  !   do m=2,maxitersSC
  !      ! calculate i (cyclical index for storing history) and pul_mx
  !      !ORI i = mod(m, maxpulaySC)
  !      !ORI if(i==0) i=maxpulaySC
  !      !ORI pul_mx = min(m, maxpulaySC)
  !      !Reset Version   TM Nov2007
  !      i = mod(m-IterPulayReset+1, maxpulaySC)
  !      if(i==0) i=maxpulaySC
  !      pul_mx = min(m-IterPulayReset+1, maxpulaySC)
  !      !if(inode==ionode) write(io_lun,fmt='(8x,"Pulay iteration ",i4)') m
  !      do j=1,n_my_grid_points
  !         rho_pul(j,i) = rho(j)
  !      end do
  !      ! Generate new charge and find residual
  !      if(iprint_SC>1) call dump_charge(rho,size,inode)
  !      call get_new_rho(.false., reset_L, fixed_potential, vary_mu, &
  !           n_L_iterations, L_tol, mu, total_energy, rho, rho1, size)
  !      !OLD if(iprint_SC>1) call dump_charge(rho,size,inode)
  !      !call dump_charge2('ch'//digitstr(m:m),rho,size,inode)
  !      do j=1,n_my_grid_points
  !         resid(j) = rho1(j) - rho(j)
  !      end do
  !      !tmp = grid_point_volume*asum(n_my_grid_points,resid,1)
  !      !write(io_lun,*) 'Sum of resid: ',tmp
  !      R0 = dot(n_my_grid_points,resid,1,resid,1)
  !      call gsum(R0)
  !      R0 = sqrt(grid_point_volume*R0)/ne_in_cell
  !      if(inode==ionode) write(io_lun,fmt='(8x,"Pulay iteration ",i5," Residual is ",e12.5,/)') m,R0
  !      !call dump_charge2('re'//digitstr(m:m),resid,size,inode)
  !      !call dump_locps(resid,size,inode)
  !      !Reset Pulay Iterations
  !      if(R0 > R0_old ) icounter_fail = icounter_fail+1
  !      if(icounter_fail > mx_fail) then
  !         if(inode == ionode) write(io_lun,*) ' Pulay iteration is reset !!  at ',m,'  th iteration'
  !         reset_Pulay = .true.
  !      endif
  !      R0_old = R0
  !      !Reset Pulay Iterations  -- introduced by TM, Nov2007
  !      do j=1,n_my_grid_points
  !         R(j,i) = resid(j)
  !      end do
  !      if(R0<self_tol) then
  !         if(inode==ionode) write(io_lun,'(8x,"Reached self-consistent tolerance")')
  !         done = .true.
  !         call dealloc_PulayMiXSCA
  !         return
  !      end if
  !      if(R0<EndLinearMixing) then
  !         if(inode==ionode) write(io_lun,'(8x,"Reached transition to LateSC")')
  !         call dealloc_PulayMixSCA
  !         return
  !      end if
  !      ! now build new rho
  !      Aij = zero
  !      do ii=1,pul_mx
  !         R1 = dot(n_my_grid_points, R(:,ii), 1, R(:,ii),1)
  !         call gsum(R1)
  !         Aij(ii,ii) = R1
  !         if(ii>1) then
  !            do j=1,ii-1
  !               R1 = dot(n_my_grid_points, R(:,ii), 1, R(:,j),1)
  !               call gsum(R1)
  !               Aij(j,ii) = R1
  !               Aij(ii,j) = R1
  !            enddo
  !         end if
  !         !if(iprint_SC>2) write(io_lun,fmt='(5f16.12)') Aij(:,ii)
  !      enddo
  !      ! Solve to get alphas
  !      call DoPulay(Aij,alph,pul_mx,maxpulaySC,inode,ionode)
  !      if(inode==ionode.AND.iprint_SC>2) write(io_lun,*) 'alph: ',alph(1:pul_mx)
  !      ! Create delta rho - add on new rho
  !      rho = 0.0_double
  !      max_neg = 0.0_double
  !      do ii=1,pul_mx
  !         Kresid = zero
  !         Kresid(1:n_my_grid_points) = R(1:n_my_grid_points,ii)
  !         call kerker(Kresid,maxngrid,q0)
  !         !if(ii==pul_mx) then
  !         rho(1:n_my_grid_points) = rho(1:n_my_grid_points) + &
  !              alph(ii)*(rho_pul(1:n_my_grid_points,ii) + A*Kresid(1:n_my_grid_points))
  !         !else
  !         ! rho(1:n_my_grid_points) = rho(1:n_my_grid_points) + &
  !         !      alph(ii)*rho_pul(1:n_my_grid_points,ii)
  !         !endif
  !      end do
  !      do j=1,n_my_grid_points
  !         if(rho(j)<0.0_double) then
  !            if(-rho(j)>max_neg) max_neg = -rho(j)
  !         end if
  !      end do
  !      if(max_neg>0.0_double.AND.iprint_SC>2.AND.inode==ionode) write(io_lun,*) i,' premix Max negative dens on node ',inode,max_neg
  !      n_iters = n_iters+1
  !      !Reset Pulay Iterations  -- introduced by TM, Nov2007
  !      if(reset_Pulay) then
  !         rho_pul(1:n_my_grid_points, 1) = rho_pul(1:n_my_grid_points, i)
  !         R(1:n_my_grid_points, 1) = R(1:n_my_grid_points, i)
  !         IterPulayReset = m
  !         icounter_fail=0
  !      endif
  !      reset_Pulay = .false.
  !      !Reset Pulay Iterations  -- introduced by TM, Nov2007
  !   end do
  !   ndone = n_iters
  !   call dealloc_PulayMixSCA
  !   return
  ! contains

  !   subroutine dealloc_PulayMixSCA

  !     implicit none

  !     integer :: stat

  !     deallocate(rho_pul, R,rho1, resid, Kresid, Aij, alph, STAT=stat)
  !     if(stat/=0) call cq_abort("Deallocation error in PulayMixSCA: ",maxngrid, maxpulaySC)
  !     call reg_dealloc_mem(area_SC, maxngrid+2*maxngrid*(1+maxpulaySC)+maxpulaySC*(maxpulaySC+1), type_dbl)
  !     return
  !   end subroutine dealloc_PulayMixSCA

  ! end subroutine PulayMixSCA
  ! !!*****


  ! !!****f* SelfCon_module/PulayMixSCB
  ! !! PURPOSE
  ! !! INPUTS
  ! !! OUTPUT
  ! !! RETURN VALUE
  ! !! AUTHOR
  ! !!   David Bowler
  ! !! CREATION DATE
  ! !!
  ! !! MODIFICATION HISTORY
  ! !!   2011/09/13 L.Tong
  ! !!     Removed obsolete dependence on number_of_bands
  ! !! SOURCE
  ! !!
  ! subroutine PulayMixSCB(done,ndone,self_tol, reset_L, &
  !      fixed_potential, vary_mu, n_L_iterations, L_tol, mu, &
  !      total_energy,rho, size)
  !   use datatypes
  !   use numbers
  !   use PosTan
  !   use Pulay, only: DoPulay
  !   use GenBlas
  !   use dimens, only: n_my_grid_points, grid_point_volume
  !   use EarlySCMod, only: get_new_rho, mixtwo
  !   use GenComms, only: gsum, cq_abort, inode, ionode
  !   use io_module, only: dump_charge, dump_charge2
  !   use hartree_module, only: kerker
  !   use global_module, only: ne_in_cell, flag_continue_on_SC_fail, flag_SCconverged
  !   use maxima_module, only: maxngrid

  !   implicit none

  !   ! Passed variables
  !   logical :: done
  !   logical :: vary_mu, fixed_potential, reset_L

  !   integer :: size
  !   integer :: ndone
  !   integer :: n_L_iterations

  !   real(double) :: self_tol
  !   real(double) :: L_tol, mixing, mu
  !   real(double) :: total_energy
  !   real(double), dimension(size) :: rho

  !   ! Local variables
  !   integer :: n_iters, n_pulay, npmod, pul_mx, i, j, ii, m, stat, thispulaySC
  !   real(double) :: R0, R1, E0, E1, dE, stLtol, tmp, max_neg

  !   real(double), allocatable, dimension(:,:) :: rho_pul, R
  !   real(double), allocatable, dimension(:) :: rho1, resid, Kresid
  !   real(double), allocatable, dimension(:,:) :: Aij
  !   real(double), allocatable, dimension(:) :: alph

  !   character(len=11) :: digitstr = "01234567890"

  !   allocate(rho_pul(maxngrid,maxpulaySC), R(maxngrid,maxpulaySC),rho1(maxngrid), &
  !        resid(maxngrid), Kresid(maxngrid), Aij(maxpulaySC, maxpulaySC), alph(maxpulaySC), STAT=stat)
  !   if(stat/=0) call cq_abort("Allocation error in PulayMixSC: ",maxngrid, maxpulaySC)
  !   rho_pul = zero
  !   R = zero
  !   rho1 = zero
  !   resid = zero
  !   Kresid = zero
  !   Aij = zero
  !   alph = zero
  !   if(inode==ionode) write(io_lun,fmt='(8x,"Starting Pulay mixing, A = ",f6.3," q0= ",f7.4)') A, q0
  !   done = .false.
  !   n_iters = ndone
  !   if(n_iters>=maxitersSC) then
  !      if(.NOT.flag_continue_on_SC_fail) call cq_abort('lateSC: too many iterations: ',n_iters, maxitersSC)
  !      flag_SCconverged = .false.
  !      done = .true.
  !      return
  !   endif
  !   ! m=1
  !   max_neg = 0.0_double
  !   do j=1,n_my_grid_points
  !      rho_pul(j,1) = rho(j)
  !   end do
  !   if(inode==ionode.AND.max_neg>0.0_double) write(io_lun,*) 'Init Max negative dens on node ',inode,max_neg
  !   call get_new_rho(.false., reset_L, fixed_potential, vary_mu, &
  !        n_L_iterations, L_tol, mu, total_energy, rho, rho1, size)
  !   ! Evaluate residual
  !   do j=1,n_my_grid_points
  !      resid(j) = rho1(j) - rho(j)
  !   end do
  !   call kerker(resid,maxngrid,q0)
  !   !tmp = grid_point_volume*rsum(n_my_grid_points,resid,1)
  !   !write(io_lun,*) 'Sum of resid: ',tmp
  !   R0 = dot(n_my_grid_points,resid,1,resid,1)
  !   call gsum(R0)
  !   R0 = sqrt(grid_point_volume*R0)/ne_in_cell
  !   if(R0<self_tol) then
  !      if(inode==ionode) write(io_lun,'(8x,"Reached tolerance")')
  !      done = .true.
  !      call dealloc_PulayMixSCB
  !      return
  !   end if
  !   if(R0<EndLinearMixing) then
  !      if(inode==ionode) write(io_lun,'(8x,"Reached transition to LateSC")')
  !      call dealloc_PulayMixSCB
  !      return
  !   end if
  !   if(inode==ionode) write(io_lun,*) 'Residual is ',R0
  !   do j=1,n_my_grid_points
  !      R(j,1) = resid(j)
  !   end do
  !   rho1 = zero
  !   rho1 = rho + A*resid
  !   !tmp = grid_point_volume*rsum(n_my_grid_points,rho1,1)
  !   !call gsum(tmp)
  !   !write(io_lun,*) 'Sum of rho1: ',tmp
  !   !rho1 = ne_in_cell*rho1/tmp
  !   rho = rho1
  !   ! Store deltarho
  !   n_iters = n_iters+1
  !   m = 1
  !   i = 1
  !   do while (m<=maxitersSC)
  !      m = m+1
  !      i = i+1
  !      ! calculate i (cyclical index for storing history) and pul_mx
  !      i = mod(i, maxpulaySC)
  !      if(i==0) i=maxpulaySC
  !      ! Need something here to replace WORST residual
  !      pul_mx = min(m, maxpulaySC)
  !      if(inode==ionode) write(io_lun,fmt='(8x,"Pulay iteration ",i4," counters are ",2i4)') m, i, pul_mx
  !      do j=1,n_my_grid_points
  !         rho_pul(j,i) = rho(j)
  !      end do
  !      ! Generate new charge and find residual
  !      call get_new_rho(.false., reset_L, fixed_potential, vary_mu, &
  !           n_L_iterations, L_tol, mu, total_energy, rho, rho1, size)
  !      !call dump_charge2('ch'//digitstr(m:m),rho,size,inode)
  !      do j=1,n_my_grid_points
  !         resid(j) = rho1(j) - rho(j)
  !      end do
  !      call kerker(resid,maxngrid,q0)
  !      R0 = dot(n_my_grid_points,resid,1,resid,1)
  !      call gsum(R0)
  !      R0 = sqrt(grid_point_volume*R0)/ne_in_cell
  !      if(inode==ionode) write(io_lun,*) 'Residual is ',R0
  !      !call dump_charge2('re'//digitstr(m:m),resid,size,inode)
  !      !call dump_locps(resid,size,inode)
  !      do j=1,n_my_grid_points
  !         R(j,i) = resid(j)
  !      end do
  !      if(R0<self_tol) then
  !         if(inode==ionode) write(io_lun,'(8x,"Reached tolerance")')
  !         done = .true.
  !         call dealloc_PulayMixSCB
  !         return
  !      end if
  !      if(R0<EndLinearMixing) then
  !         if(inode==ionode) write(io_lun,'(8x,"Reached transition to LateSC")')
  !         call dealloc_PulayMixSCB
  !         return
  !      end if
  !      ! now build new rho
  !      Aij = zero
  !      if(inode==ionode) write(io_lun,*) 'A is :'
  !      do ii=1,pul_mx
  !         do j=1,pul_mx
  !            R1 = dot(n_my_grid_points, R(:,ii), 1, R(:,j),1)
  !            call gsum(R1)
  !            Aij(j,ii) = R1
  !         enddo
  !         if(inode==ionode) write(io_lun,fmt='(5f16.12)') Aij(:,ii)
  !      enddo
  !      !if(minA/maxA<0.001_double) then
  !      !   write(io_lun,*) '
  !      !end if
  !      ! Solve to get alphas
  !      call DoPulay(Aij,alph,pul_mx,maxpulaySC,inode,ionode)
  !      if(inode==ionode)write(io_lun,*) 'alph: ',alph(1:pul_mx)
  !      ! Create delta rho - add on new rho
  !      rho = 0.0_double
  !      max_neg = 0.0_double
  !      do ii=1,pul_mx
  !         !Kresid = zero
  !         !Kresid(1:n_my_grid_points) = R(1:n_my_grid_points,ii)
  !         !call kerker(Kresid,maxngrid,q0)
  !         rho(1:n_my_grid_points) = rho(1:n_my_grid_points) + &
  !              alph(ii)*(rho_pul(1:n_my_grid_points,ii) + A*R(1:n_my_grid_points,ii))
  !         !alph(ii)*(rho_pul(1:n_my_grid_points,ii) + A*Kresid(1:n_my_grid_points))
  !      end do
  !      tmp = grid_point_volume*rsum(n_my_grid_points,rho,1)
  !      call gsum(tmp)
  !      if(inode==ionode)write(io_lun,*) 'Sum of rho1: ',tmp
  !      rho = ne_in_cell*rho/tmp
  !      !rho = rho*ne_in_cell/tmp
  !      !do j=1,n_my_grid_points
  !      !if(rho(j)<0.0_double) then
  !      !   if(-rho(j)>max_neg) max_neg = -rho(j)
  !      !   rho(j) = 0.0_double
  !      !end if
  !      !end do
  !      if(max_neg>0.0_double.AND.inode==ionode) write(io_lun,*) i,' premix Max negative dens on node ',inode,max_neg
  !      n_iters = n_iters+1
  !   end do
  !   ndone = n_iters
  !   call dealloc_PulayMixSCB
  !   return
  ! contains

  !   subroutine dealloc_PulayMixSCB

  !     use GenComms, only: cq_abort
  !     implicit none
  !     integer :: stat

  !     deallocate(rho_pul, R,rho1, resid, Kresid, Aij, alph, STAT=stat)
  !     if(stat/=0) call cq_abort("Deallocation error in PulayMixSCB: ",maxngrid, maxpulaySC)
  !     return
  !   end subroutine dealloc_PulayMixSCB

  ! end subroutine PulayMixSCB
  ! !!*****


  ! !!****f* SelfCon_module/PulayMixSCC *
  ! !!
  ! !!  NAME
  ! !!   PulayMixSCC -- PulayMixSC version C, Lianheng's version
  ! !!  USAGE
  ! !!   CALL PulayMixSCC (done,ndone,self_tol, reset_L, fixed_potential, &
  ! !!                     vary_mu, n_L_iterations, number_of_bands, &
  ! !!                     L_tol, mu, total_energy, rho, size)
  ! !!  PURPOSE
  ! !!   Performs Pulay mixing of charge density with Kerker
  ! !!   preconditioning and wave dependent metric in Conquest
  ! !!  INPUTS
  ! !!
  ! !!  USES
  ! !!
  ! !!  AUTHOR
  ! !!   Lianheng Tong
  ! !!  CREATION DATE
  ! !!   2010/07/27
  ! !!  MODIFICATION HISTORY
  ! !!   2011/09/13 L.Tong
  ! !!     Removed obsolete dependence on number_of_bands
  ! !!  SOURCE
  ! !!
  ! subroutine PulayMixSCC (done,ndone,self_tol, reset_L, &
  !      fixed_potential, vary_mu, n_L_iterations, L_tol, mu, &
  !      total_energy,rho, size)
  !   use datatypes
  !   use numbers
  !   use PosTan
  !   use Pulay, only: DoPulay
  !   use GenBlas
  !   use dimens, only: n_my_grid_points, grid_point_volume
  !   use EarlySCMod, only: get_new_rho, mixtwo
  !   use GenComms, only: gsum, cq_abort, inode, ionode
  !   use io_module, only: dump_charge, dump_charge2
  !   use hartree_module, only: kerker, kerker_and_wdmetric, wdmetric
  !   use global_module, only: ne_in_cell, iprint_SC, area_SC, &
  !        flag_continue_on_SC_fail, flag_SCconverged
  !   use memory_module, only: reg_alloc_mem, reg_dealloc_mem, type_dbl
  !   use maxima_module, only: maxngrid

  !   implicit none

  !   ! Passed variables
  !   logical :: done, vary_mu, fixed_potential, reset_L
  !   integer :: size, ndone, n_L_iterations
  !   real(double) :: self_tol, L_tol, mixing, mu, total_energy
  !   real(double), dimension(size) :: rho

  !   ! Local variables
  !   integer :: n_iters, pul_mx, i, j, ii, m, stat
  !   real(double) :: R0, R1
  !   real(double), allocatable, dimension(:,:) :: rho_pul, R_pul
  !   real(double), allocatable, dimension(:) :: rho_out, resid
  !   real(double), allocatable, dimension(:,:) :: Aij
  !   real(double), allocatable, dimension(:) :: alpha
  !   ! Kerker only
  !   real(double), allocatable, dimension(:,:) :: KR_pul
  !   ! wdmetric only
  !   real(double), allocatable, dimension(:) :: resid_cov
  !   real(double), allocatable, dimension(:,:) :: Rcov_pul

  !   !Reset Pulay Iterations  -- introduced by TM, Nov2007
  !   integer, parameter :: mx_fail = 3
  !   integer :: IterPulayReset=1, icounter_fail=0
  !   logical :: reset_Pulay=.false.
  !   real(double) :: R0_old

  !   allocate (rho_pul(maxngrid,maxpulaySC), R_pul(maxngrid,maxpulaySC), rho_out(maxngrid), &
  !        resid(maxngrid), Aij(maxpulaySC, maxpulaySC), alpha(maxpulaySC), STAT=stat)
  !   if (stat /= 0) call cq_abort ("SelfCon/PulayMixSCC: Allocation error: ", maxngrid, maxpulaySC)
  !   call reg_alloc_mem (area_SC, 2*maxngrid*maxpulaySC+2*maxngrid+maxpulaySC*(maxpulaySC+1), type_dbl)
  !   rho_pul = zero
  !   R_pul = zero
  !   rho_out = zero
  !   resid = zero
  !   Aij = zero
  !   alpha = zero
  !   ! for Kerker preconditioning
  !   if (flag_Kerker) then
  !      allocate (KR_pul(maxngrid,maxpulaySC), STAT=stat)
  !      if (stat /= 0) call cq_abort ("SelfCon/PulayMixSCC: Allocation error, KR_pul: ", maxngrid, maxpulaySC)
  !      call reg_alloc_mem (area_SC, maxngrid*maxpulaySC, type_dbl)
  !      KR_pul = zero
  !   end if
  !   ! for wave dependent metric method
  !   if (flag_wdmetric) then
  !      allocate (Rcov_pul(maxngrid,maxpulaySC), resid_cov(maxngrid), STAT=stat)
  !      if (stat /= 0) call cq_abort ("SelfCon/PulayMixSCC: Allocation error, Rcov_pul, resid_cov: ",maxngrid, maxpulaySC)
  !      call reg_alloc_mem (area_SC, maxngrid*(maxpulaySC+1), type_dbl)
  !      Rcov_pul = zero
  !      resid_cov = zero
  !   end if

  !   ! write out start information
  !   if (inode == ionode) then
  !      write (io_lun, fmt='(8x,"Starting Pulay mixing, A = ",f6.3)') A
  !      if (flag_Kerker) write (io_lun, fmt='(10x,"with Kerker preconditioning, q0 = ", f6.3)') q0
  !      if (flag_wdmetric) write (io_lun, fmt='(10x,"with wave dependent metric, q1 = ", f6.3)') q1
  !   end if

  !   done = .false.
  !   n_iters = ndone
  !   m = 1
  !   if (n_iters >= maxitersSC) then
  !      if (.NOT.flag_continue_on_SC_fail) call cq_abort ('SelfCon/PulayMixSCC: too many SCF iterations: ', n_iters, maxitersSC)
  !      flag_SCconverged = .false.
  !      done = .true.
  !      return
  !   end if

  !   ! store the first rho_in in pulay history
  !   do j = 1, n_my_grid_points
  !      rho_pul(j,1) = rho(j)
  !   end do
  !   ! from rho_in, get rho_out
  !   call get_new_rho (.false., reset_L, fixed_potential, vary_mu, &
  !        n_L_iterations, L_tol, mu, total_energy, rho, rho_out, size)
  !   ! Evaluate residue resid = rho_out - rho_in
  !   do j = 1, n_my_grid_points
  !      resid(j) = rho_out(j) - rho(j)
  !   end do
  !   ! calculate the norm of residue
  !   R0 = dot (n_my_grid_points, resid, 1, resid, 1)
  !   ! gather for all processor nodes
  !   call gsum (R0)
  !   ! normalise
  !   R0 = sqrt (grid_point_volume*R0)/ne_in_cell
  !   if (R0 < self_tol) then
  !      if (inode == ionode) write (io_lun, fmt='(8x,"Reached self-&
  !           &consistency tolerance")')
  !      done = .true.
  !      call dealloc_PulayMiXSCC
  !      return
  !   end if
  !   if (R0 < EndLinearMixing) then
  !      if (inode == ionode) write (io_lun, fmt='(8x,"Reached &
  !           &transition to LateSC")')
  !      call dealloc_PulayMiXSCC
  !      return
  !   end if
  !   if (inode == ionode) write (io_lun, fmt='(8x,"Pulay iteration ",&
  !        &i5," Residual is ",e12.5,/)') m, R0
  !   ! store resid to pulay history
  !   do j = 1, n_my_grid_points
  !      R_pul(j,1) = resid(j)
  !   end do
  !   ! do Kerker preconditioning if required
  !   if (flag_Kerker) then
  !      ! get wave dependent metric convariant R as well if required
  !      if (flag_wdmetric) then
  !         call kerker_and_wdmetric (resid, resid_cov, maxngrid, q0, q1)
  !         ! store the preconditioned and covariant version of residue
  !         ! to histiry
  !         do  j = 1, n_my_grid_points
  !            KR_pul(j,1) = resid(j)
  !            Rcov_pul(j,1) = resid_cov(j)
  !         end do
  !      else
  !         call kerker (resid, maxngrid, q0)
  !         ! store the precondiioned residue to history
  !         do j = 1, n_my_grid_points
  !            KR_pul(j,1) = resid(j)
  !         end do
  !      end if
  !   else
  !      ! do wave dependent metric without kerker preconditioning if
  !      ! required
  !      if (flag_wdmetric) then
  !         call wdmetric (resid, resid_cov, maxngrid, q1)
  !         ! store the covariant version of residue
  !         do j = 1, n_my_grid_points
  !            Rcov_pul(j,1) = resid_cov(j)
  !         end do
  !      end if
  !   end if
  !   ! Do mixing
  !   call axpy (n_my_grid_points, A, resid, 1, rho, 1)

  !   ! finished the first SCF iteration
  !   n_iters = n_iters + 1
  !   ! Routines for Resetting Pulay Iterations (TM)
  !   R0_old = R0
  !   IterPulayReset = 1
  !   icounter_fail = 0
  !   ! Start SCF loop from second iteration
  !   SCF: do m = 2, maxitersSC
  !      ! calculate i (cyclical index for storing history), including TM's Pulay reset
  !      i = mod (m-IterPulayReset+1, maxpulaySC)
  !      if (i == 0) i = maxpulaySC
  !      ! calculated the number of pulay histories stored
  !      pul_mx = min (m-IterPulayReset+1, maxpulaySC)
  !      ! store the updated density from the previous step to history
  !      do j = 1, n_my_grid_points
  !         rho_pul(j,i) = rho(j)
  !      end do
  !      ! print out charge
  !      if (iprint_SC>1) call dump_charge (rho,size,inode)
  !      ! genertate new charge density
  !      call get_new_rho (.false., reset_L, fixed_potential, vary_mu, &
  !           n_L_iterations, L_tol, mu, total_energy, rho, rho_out, &
  !           size)
  !      ! get residue
  !      do j = 1, n_my_grid_points
  !         resid(j) = rho_out(j) - rho(j)
  !      end do
  !      R0 = dot (n_my_grid_points, resid, 1, resid, 1)
  !      call gsum(R0)
  !      R0 = sqrt (grid_point_volume*R0)/ne_in_cell
  !      ! print out SC iteration information
  !      if (inode == ionode) write (io_lun, fmt='(8x,"Pulay iteration &
  !           &",i5," Residual is ",e12.5,/)') m, R0
  !      ! check if Pulay SC has converged
  !      if (R0 < self_tol) then
  !         if (inode == ionode) write (io_lun, fmt='(8x,"Reached self-&
  !              &consistent tolerance")')
  !         done = .true.
  !         call dealloc_PulayMiXSCC
  !         return
  !      end if
  !      if (R0 < EndLinearMixing) then
  !         if (inode == ionode) write (io_lun, fmt='(8x,"Reached &
  !              &transition to LateSC")')
  !         call dealloc_PulayMixSCC
  !         return
  !      end if
  !      ! check if Pulay SC needs to be reset
  !      if (R0 > R0_old) icounter_fail = icounter_fail+1
  !      if (icounter_fail > mx_fail) then
  !         if (inode == ionode) write(io_lun, *) ' Pulay iteration is  &
  !              &              reset !!  at ', m, '  th iteration'
  !         reset_Pulay = .true.
  !      endif
  !      R0_old = R0
  !      ! continue on with SC, store resid to history
  !      ! but before doing so, remember the old R_pul(j,i) to takeaway
  !      ! from Kerker_sum, Kerker only
  !      do j = 1, n_my_grid_points
  !         R_pul(j,i) = resid(j)
  !      end do
  !      ! calculate new rho
  !      ! with kerker preconditioning
  !      if (flag_Kerker) then
  !         ! with wave dependent metric
  !         if (flag_wdmetric) then
  !            call kerker_and_wdmetric (resid, resid_cov, maxngrid, q0,&
  !                 q1)
  !            ! store preconditioned and covariant residue to history
  !            do j = 1, n_my_grid_points
  !               KR_pul(j,i) = resid(j)
  !               Rcov_pul(j,i) = resid_cov(j)
  !            end do
  !            ! get alpha_i
  !            Aij = zero
  !            do ii = 1, pul_mx
  !               ! diagonal elements of Aij
  !               R1 = dot (n_my_grid_points, Rcov_pul(:,ii), 1, &
  !                    R_pul(:,ii), 1)
  !               call gsum (R1)
  !               Aij(ii,ii) = R1
  !               ! the rest
  !               if (ii > 1) then
  !                  do j = 1, ii-1
  !                     R1 = dot (n_my_grid_points, Rcov_pul(:,ii), 1, &
  !                          R_pul(:,j), 1)
  !                     call gsum (R1)
  !                     Aij(ii,j) = R1
  !                     ! Aij is symmetric even with wave dependent
  !                     ! metric
  !                     Aij(j,ii) = R1
  !                  end do
  !               end if
  !            end do
  !         else  ! If no wave dependent metric
  !            call kerker (resid, maxngrid, q0)
  !            ! store preconditioned residue to history
  !            do j = 1, n_my_grid_points
  !               KR_pul(j,i) = resid(j)
  !            end do
  !            ! get alpha_i
  !            Aij = zero
  !            do ii = 1, pul_mx
  !               ! diagonal elements of Aij
  !               R1 = dot (n_my_grid_points, R_pul(:,ii), 1, R_pul(:,&
  !                    ii), 1)
  !               call gsum (R1)
  !               Aij(ii,ii) = R1
  !               ! the rest
  !               if (ii > 1) then
  !                  do j = 1, ii-1
  !                     R1 = dot (n_my_grid_points, R_pul(:,ii), 1, &
  !                          R_pul(:,j), 1)
  !                     call gsum (R1)
  !                     Aij(ii,j) = R1
  !                     Aij(j,ii) = R1
  !                  end do
  !               end if
  !            end do
  !         end if ! if (flag_wdmetric)
  !         ! solve alpha(i) = sum_j Aji^-1 / sum_ij Aji^-1
  !         call DoPulay (Aij, alpha, pul_mx, maxpulaySC, inode, &
  !              ionode)
  !         ! Do mixing
  !         rho = zero
  !         do ii = 1, pul_mx
  !            call axpy (n_my_grid_points, alpha(ii), rho_pul(:,ii), 1,&
  !                 rho, 1)
  !            call axpy (n_my_grid_points, alpha(ii)*A, KR_pul(:,ii), &
  !                 1, rho, 1)
  !         end do
  !      else  ! if no Kerker
  !         ! do wave dependent metric without kerker preconditioning
  !         ! if required
  !         if (flag_wdmetric) then
  !            call wdmetric (resid, resid_cov, maxngrid, q1)
  !            ! store the covariant version of residue to history
  !            do j = 1, n_my_grid_points
  !               Rcov_pul(j,i) = resid_cov(j)
  !            end do
  !            ! get alpha_i
  !            Aij = zero
  !            do ii = 1, pul_mx
  !               ! diagonal elements of Aij
  !               R1 = dot (n_my_grid_points, Rcov_pul(:,ii), 1, &
  !                    R_pul(:,ii), 1)
  !               call gsum (R1)
  !               Aij(ii,ii) = R1
  !               ! the rest
  !               if (ii > 1) then
  !                  do j = 1, ii-1
  !                     R1 = dot (n_my_grid_points, Rcov_pul(:,ii), 1, &
  !                          R_pul(:,j), 1)
  !                     call gsum (R1)
  !                     Aij(ii,j) = R1
  !                     ! Aij is symmetric even with wave dependent
  !                     ! metric
  !                     Aij(j,ii) = R1
  !                  end do
  !               end if
  !            end do
  !         else ! if no Kerker and no wave dependent metric
  !            ! get alpha_i
  !            Aij = zero
  !            do ii = 1, pul_mx
  !               ! diagonal elements of Aij
  !               R1 = dot (n_my_grid_points, R_pul(:,ii), 1, R_pul(:,&
  !                    ii), 1)
  !               call gsum (R1)
  !               Aij(ii,ii) = R1
  !               ! the rest
  !               if (ii > 1) then
  !                  do j = 1, ii-1
  !                     R1 = dot (n_my_grid_points, R_pul(:,ii), 1, &
  !                          R_pul(:,j), 1)
  !                     call gsum (R1)
  !                     Aij(ii,j) = R1
  !                     ! Aij is symmetric even with wave dependent
  !                     ! metric
  !                     Aij(j,ii) = R1
  !                  end do
  !               end if
  !            end do
  !         end if ! if (flag_wdmetric)
  !         ! solve alpha(i) = sum_j Aji^-1 / sum_ij Aji^-1
  !         call DoPulay (Aij, alpha, pul_mx, maxpulaySC, inode, &
  !              ionode)
  !         ! Do mixing
  !         rho = zero
  !         do ii = 1, pul_mx
  !            call axpy (n_my_grid_points, alpha(ii), rho_pul(:,ii), 1,&
  !                 rho, 1)
  !            call axpy (n_my_grid_points, alpha(ii)*A, R_pul(:,ii), 1,&
  !                 rho, 1)
  !         end do
  !      end if ! if (flag_Kerker)
  !      n_iters = n_iters + 1
  !      ! Reset Pulay iterations if required
  !      if (reset_Pulay) then
  !         rho_pul(1:n_my_grid_points, 1) = rho_pul(1:n_my_grid_points,&
  !              i)
  !         R_pul(1:n_my_grid_points, 1) = R_pul(1:n_my_grid_points, i)
  !         if (flag_Kerker) KR_pul(1:n_my_grid_points, 1) = &
  !              KR_pul(1:n_my_grid_points, i)
  !         if (flag_wdmetric) Rcov_pul(1:n_my_grid_points, 1) = &
  !              Rcov_pul(1:n_my_grid_points, i)
  !         IterPulayReset = m
  !         icounter_fail = 0
  !      end if
  !      reset_Pulay = .false.
  !   end do SCF
  !   ndone = n_iters
  !   call dealloc_PulayMixSCC
  !   return

  ! contains

  !   subroutine dealloc_PulayMixSCC

  !     implicit none
  !     integer :: stat

  !     deallocate (rho_pul, R_pul, rho_out, resid, Aij, alpha, STAT=&
  !          stat)
  !     if (stat /= 0) call cq_abort ("Deallocation error in &
  !          &PulayMixSCC: ", maxngrid, maxpulaySC)
  !     call reg_dealloc_mem (area_SC, 2*maxngrid*maxpulaySC+2*maxngrid+&
  !          maxpulaySC*(maxpulaySC+1), type_dbl)

  !     if (flag_Kerker) then
  !        deallocate (KR_pul, STAT=stat)
  !        if (stat /= 0) call cq_abort ("SelfCon/dealloc_PulayMixSCC: &
  !             &Deallocation error, KR_pul: ", maxngrid, maxpulaySC)
  !        call reg_dealloc_mem (area_SC, maxngrid*maxpulaySC, type_dbl)
  !     end if

  !     if (flag_wdmetric) then
  !        deallocate (Rcov_pul, resid_cov, STAT=stat)
  !        if (stat /= 0) call cq_abort ("SelfCon/dealloc_PulayMixSCC: &
  !             &Deallocation error, Rcov_pul, resid_cov: ", maxngrid, &
  !             maxpulaySC)
  !        call reg_dealloc_mem (area_SC, maxngrid*(maxpulaySC+1), &
  !             type_dbl)
  !     end if

  !     return

  !   end subroutine dealloc_PulayMixSCC

  ! end subroutine PulayMixSCC


  !!****f* SelfCon_module/PulayMixSC_spin
  !! PURPOSE
  !!   Performs Pulay mixing of charge density with Kerker
  !!   preconditioning and wavedependent metric in Conquest for spin
  !!   polarised calculations.
  !!
  !!   Based on SelfCon_module/PulayMixSCC
  !!
  !!   Note that the Pulay coefficients alpha(i) and the metric Aij
  !!   should be calculated from the total Residual R = R^up + R^dn and
  !!   the total density rho = rho^up + rho^dn. This is because for
  !!   Pulay method to work the alpha(i) are constrained so that
  !!
  !!             \sum_i \alpha(i) = 1
  !!
  !!   This conserves the electron numbers IF the input densities
  !!   them-selves all conserve electron number. However for the case
  !!   where the spin populations are allowed to change, rho^up and
  !!   rho^dn no longer conserve electron numbers. Hence if we have
  !!   separate alpha's for different spin channels. The constriant
  !!
  !!            \sum_i \alpha^sigma(i) = 1
  !!
  !!   is no longer sufficient for fixing electron numbers in each spin
  !!   channels and the total electron number may vary as a result.
  !!
  !! INPUTS
  !! OUTPUT
  !! RETURN VALUE
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2011/08/30
  !! MODIFICATION HISTORY
  !!  2011/09/13 L.Tong
  !!    removed obsolete dependence on number_of_bands
  !!  2011/09/18 L.Tong
  !!    Major correction, now Pulay coefficients are spin dependent, and
  !!    different constraints for fixed or non-fixed spin populations
  !!  2012/03/20 L.Tong
  !!  - Major change to spin implementation
  !!  - Removed redundant input parameter real(double) mu
  !!  - Made R0 to be the total magnitude including both spin
  !!    channels (note the cross term)
  !!  2012/05/25 L.Tong
  !!  - Major rewrite of the subroutine
  !!  - The pulay mixing step (calculation of alpha) are now put into
  !!    a separate subroutine
  !!   2014/09/28 L.Truflandier
  !!  - Added exx_pulay_r0 to control EXX accuracy during the SCF based on
  !!    Pulay mixing: the point is to extract R0 from PulayMixSC and use it 
  !!    in get_X_matrix
  !! SOURCE
  !!
  subroutine PulayMixSC_spin(done, ndone, self_tol, reset_L, &
                             fixed_potential, vary_mu, n_L_iterations,&
                             L_tol, total_energy, rho, size)
    use datatypes
    use numbers
    use PosTan
    use Pulay,          only: DoPulay
    use GenBlas
    use dimens,         only: n_my_grid_points, grid_point_volume
    use EarlySCMod,     only: get_new_rho
    use GenComms,       only: gsum, cq_abort, inode, ionode, my_barrier
    use io_module,      only: dump_charge
    use hartree_module, only: kerker, kerker_and_wdmetric, wdmetric
    use global_module,  only: ne_in_cell, iprint_SC, area_SC,  &
                              flag_continue_on_SC_fail,        &
                              flag_SCconverged,                &
                              flag_fix_spin_population, nspin, &
                              spin_factor
    use memory_module,  only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use maxima_module,  only: maxngrid

    implicit none

    ! Passed variables
    logical      :: done, vary_mu, fixed_potential, reset_L
    integer      :: size, ndone, n_L_iterations
    real(double) :: self_tol, L_tol, mixing, total_energy
    real(double), dimension(:,:) :: rho

    ! Local variables
    integer      :: n_iters, pul_mx, ii, iter, iPulay, stat
    real(double) :: R0
    real(double), dimension(:),     allocatable :: rho_tot
    real(double), dimension(:,:,:), allocatable :: rho_pul
    real(double), dimension(:,:,:), allocatable, target :: R_pul
    ! Kerker only
    real(double), pointer, dimension(:,:,:) :: KR_pul
    ! wdmetric only
    real(double), pointer, dimension(:,:,:) :: Rcov_pul

    !Reset Pulay Iterations  -- introduced by TM, Nov2007
    integer, parameter :: mx_fail  = 3
    integer      :: IterPulayReset = 1
    integer      :: icounter_fail  = 0
    logical      :: reset_Pulay    = .false.
    integer      :: spin
    real(double) :: R0_old

    type(cq_timer)    :: tmr_std_loc

!****lat<$
    call start_timer(t=tmr_std_loc,who='PulayMixSC_spin',&
         where=5,level=1,echo=.true.)
!****lat>$

    ! allocate memories
    call allocate_PulayMixSC_spin

    ! initialise
    rho_pul = zero
    R_pul   = zero

    ! write out start information
    if (inode == ionode) then
       write (io_lun, '(8x,a,f6.3,a,f6.3)') &
             'Starting Pulay mixing, A_up = ', A(1), ' A_dn = ', A(2)
       if (nspin == 2) then
          if (flag_fix_spin_population) then
             write (io_lun, '(8x,"Spin populations are fixed.")')
          else
             write (io_lun, '(8x,"Spin populations are to be relaxed.")')
          end if
       else
          write (io_lun, '(8x,"Spin non-polarised calculation.")')
       end if
       if (flag_Kerker) &
            write (io_lun, '(10x,a,f6.3)') &
                  'with Kerker preconditioning, q0 = ', q0
       if (flag_wdmetric) &
            write (io_lun, '(10x,a,f6.3)') &
                  'with wave dependent metric, q1 = ', q1
    end if

    ! set counters
    done    = .false.
    n_iters = ndone
    iter    = 1
    if (n_iters >= maxitersSC) then
       if (.not. flag_continue_on_SC_fail) &
            call cq_abort('SelfCon/PulayMixSC_spin: too many SCF iterations: ', &
                          n_iters, maxitersSC)
       flag_SCconverged = .false.
       done = .true.
       return
    end if

    ! store rho and calculate residuals and store in pulay history slot 1
    call update_pulay_history(1, rho, reset_L, fixed_potential,     &
                              vary_mu, n_L_iterations, L_tol,       &
                              total_energy, rho_pul, R_pul, KR_pul, &
                              Rcov_pul)

    ! Evaluate magnitute of residual, note do not include cross terms
    R0 = zero
    do spin = 1, nspin
       R0 = R0 + spin_factor * &
            dot(n_my_grid_points, R_pul(:,1,spin), 1, R_pul(:,1,spin), 1)
    end do
    ! cross term
    ! R0 = R0 + two * &
    !      dot(n_my_grid_points, R_pul(:,1,1), 1, R_pul(:,1,nspin), 1)
    call gsum(R0)
    R0 = sqrt(grid_point_volume * R0) / ne_in_cell

    ! print residual information
    if (inode == ionode) then
       write (io_lun, '(8x,a,i5,a,e12.5)') &
             'Pulay iteration ', iter, ' Residual is ', R0
    end if

    ! check if they have reached tolerance
    if (R0 < self_tol) then
       if (inode == ionode) &
            write (io_lun, '(8x,"Reached self-consistency tolerance")')
       done = .true.
       call deallocate_PulayMiXSC_spin
       return
    end if
    if (R0 < EndLinearMixing) then
       if (inode == ionode) &
            write (io_lun, '(8x,"Reached transition to LateSC")')
       call deallocate_PulayMiXSC_spin
       return
    end if

    ! Do linear mixing
    do spin = 1, nspin
       call axpy(n_my_grid_points, A(spin), KR_pul(:,1,spin), 1, &
                 rho(:,spin), 1)
    end do

    ! finished the first SCF iteration (no need to optimise)
    n_iters = n_iters + 1

    ! Record Pulay Resetting information
    R0_old = R0
    IterPulayReset = 1
    icounter_fail = 0

    ! do SCF loop
    do iter = 2, maxitersSC

       ! calculate cyclic index for storing pulay history
       iPulay = mod(iter - IterPulayReset + 1, maxpulaySC)
       if (iPulay == 0) iPulay = maxpulaySC
       ! calculated the number of pulay histories stored
       pul_mx = min(iter - IterPulayReset + 1, maxpulaySC)

       ! print out charge
       if (iprint_SC > 1) then
          if (nspin == 1) then
             rho_tot(:) = spin_factor * rho(:,1)
             call dump_charge(rho_tot, n_my_grid_points, inode, spin=0)
          else
             call dump_charge(rho(:,1), n_my_grid_points, inode, spin=1)
             call dump_charge(rho(:,2), n_my_grid_points, inode, spin=2)
          end if
       end if

       ! Calcualate residuals and update pulay history (store in iPulay-th slot)
       call update_pulay_history(iPulay, rho, reset_L, fixed_potential, &
                                 vary_mu, n_L_iterations, L_tol,        &
                                 total_energy, rho_pul, R_pul, KR_pul,  &
                                 Rcov_pul)

       ! Evaluate magnitute of residual, note no cross terms
       R0 = zero
       do spin = 1, nspin
          R0 = R0 + spin_factor * &
               dot(n_my_grid_points, R_pul(:,iPulay,spin), 1, &
                   R_pul(:,iPulay,spin), 1)
       end do
       ! cross term
       ! R0 = R0 + two * &
       !      dot(n_my_grid_points, R_pul(:,1,1), 1, R_pul(:,1,nspin), 1)
       call gsum(R0)
       R0 = sqrt(grid_point_volume * R0) / ne_in_cell

       ! print residual information
       if (inode == ionode) then
          write (io_lun, '(8x,a,i5,a,e12.5)') &
               'Pulay iteration ', iter, ' Residual is ', R0
       end if

       ! check if they have reached tolerance
       if (R0 < self_tol) then
          if (inode == ionode) &
               write (io_lun, '(8x,"Reached self-consistency tolerance")')
          done = .true.
          call deallocate_PulayMiXSC_spin
          return
       end if
       if (R0 < EndLinearMixing) then
          if (inode == ionode) &
               write (io_lun, '(8x,"Reached transition to LateSC")')
          call deallocate_PulayMiXSC_spin
          return
       end if

       ! check if Pulay SC needs to be reset
       if (R0 > R0_old) &
            icounter_fail = icounter_fail + 1
       if (icounter_fail > mx_fail) then
          if (inode == ionode) &
               write (io_lun, *) ' Pulay iteration is reset !!  at ', iter, &
                                 ' th iteration'
          reset_Pulay = .true.
       end if
       R0_old = R0

       ! get optimum rho mixed from rho_pul, R_pul and Rcov_pul
       call get_pulay_optimal_rho(iPulay, rho, pul_mx, A, rho_pul, &
                                  R_pul, KR_pul, Rcov_pul)

       ! increment iteration counter
       n_iters = n_iters + 1

       ! Reset Pulay iterations if required
       if (reset_Pulay) then
          do spin = 1, nspin
             rho_pul(1:n_my_grid_points,1,spin) =             &
                  rho_pul(1:n_my_grid_points,iPulay,spin)
             R_pul(1:n_my_grid_points,1,spin) =               &
                  R_pul(1:n_my_grid_points,iPulay,spin)
             if (flag_Kerker) then
                KR_pul(1:n_my_grid_points,1,spin) =           &
                     KR_pul(1:n_my_grid_points,iPulay,spin)
             end if
             if (flag_wdmetric) then
                Rcov_pul(1:n_my_grid_points,1,spin) =         &
                     Rcov_pul(1:n_my_grid_points,iPulay,spin)
             end if
          end do ! spin
          IterPulayReset = iter
          icounter_fail  = 0
          reset_Pulay    = .false.
       end if

    end do ! iter (SCF)

    ndone = n_iters
    call deallocate_PulayMixSC_spin

!****lat<$
    call stop_timer(t=tmr_std_loc,who='PulayMixSC_spin',echo=.true.)
!****lat>$

    return
  contains

    subroutine allocate_PulayMixSC_spin
      implicit none
      integer :: stat
      allocate(rho_tot(maxngrid), rho_pul(maxngrid,maxpulaySC,nspin), &
               R_pul(maxngrid,maxpulaySC,nspin), STAT=stat)
      if (stat /= 0) &
           call cq_abort("PulayMixSC_spin: Error alloc mem: ", &
                         maxngrid, maxpulaySC)
      call reg_alloc_mem(area_SC, (2*maxpulaySC*nspin+1)*maxngrid, type_dbl)
      ! for Kerker preconditioning
      if (flag_Kerker) then
         allocate(KR_pul(maxngrid,maxpulaySC,nspin), STAT=stat)
         if (stat /= 0) &
              call cq_abort("SelfCon/PulayMixSCC: Allocation error, KR_pul: ", &
                            stat)
         call reg_alloc_mem(area_SC, nspin * maxngrid * maxpulaySC, type_dbl)
         KR_pul = zero
      else
         KR_pul => R_pul
      end if
      ! for wave dependent metric method
      if (flag_wdmetric) then
         allocate(Rcov_pul(maxngrid,maxpulaySC,nspin), STAT=stat)
         if (stat /= 0) &
              call cq_abort("SelfCon/PulayMixSCC: Allocation error, &
                            &Rcov_pul: ", stat)
         call reg_alloc_mem(area_SC, nspin * maxngrid * maxpulaySC, &
                            type_dbl)
         Rcov_pul = zero
      else
         Rcov_pul => R_pul
      end if
    end subroutine allocate_PulayMixSC_spin

    subroutine deallocate_PulayMixSC_spin
      implicit none
      integer :: stat
      ! For Kerker preconditioning
      if (flag_Kerker) then
         deallocate(KR_pul, STAT=stat)
         if (stat /= 0) &
              call cq_abort("dealloc_PulayMixSC_spin: Deallocation error, &
                             &KR_pul ", stat)
         call reg_dealloc_mem(area_SC, nspin * maxngrid * maxpulaySC, type_dbl)
      end if
      ! For wave dependent metric
      if (flag_wdmetric) then
         deallocate(Rcov_pul, STAT=stat)
         if (stat /= 0) &
              call cq_abort ("dealloc_PulayMixSC_spin:&
                              & Deallocation error, Rcov_pul", stat)
         call reg_dealloc_mem(area_SC, nspin * maxngrid * (maxpulaySC + 1), &
                              type_dbl)
      end if
      deallocate(rho_tot, rho_pul, R_pul, STAT=stat)
      if (stat /= 0) call cq_abort("PulayMixSC_spin: Error dealloc mem")
      call reg_dealloc_mem(area_SC, (2*maxpulaySC*nspin+1)*maxngrid, type_dbl)
    end subroutine deallocate_PulayMixSC_spin

  end subroutine PulayMixSC_spin
  !*****


  !!****f* SelfCon_module/get_pulay_optimal
  !! PURPOSE
  !!   From a given pulay history calculated optimised and mixed
  !!   density
  !! USAGE
  !! INPUTS
  !! OUTPUT
  !! RETURN VALUE
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/05/21
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  subroutine get_pulay_optimal_rho(iPulay, rho_opt, pul_max,       &
                                   mix_factor, rho_pul, resid_pul, &
                                   k_resid_pul, cov_resid_pul)
    use datatypes
    use GenBlas
    use global_module,  only: nspin, flag_fix_spin_population
    use GenComms,       only: gsum, inode, ionode
    use Pulay,          only: DoPulay
    use dimens,         only: n_my_grid_points
    use maxima_module,  only: maxngrid

    implicit none
    ! Passed parameters
    integer, intent(in) :: pul_max, iPulay
    real(double), dimension(nspin),          intent(in)  :: mix_factor
    real(double), dimension(maxngrid,nspin), intent(out) :: rho_opt
    real(double), dimension(maxngrid,maxpulaySC,nspin),  intent(in)  :: &
         rho_pul, resid_pul, k_resid_pul, cov_resid_pul
    ! Local variables
    integer      :: spin, ii, jj
    real(double) :: ne, RR
    real(double), dimension(maxpulaySC,maxpulaySC,nspin) :: Aij
    real(double), dimension(maxpulaySC,nspin) :: alpha
    type(cq_timer) :: tmr_std_loc

!****lat<$
    call start_timer(t=tmr_std_loc,who='get_pulay_optimal_rho',&
         where=5,level=2,echo=.true.)
!****lat>$

    ! calculated Aij
    Aij = zero
    do spin = 1, nspin
       do ii = 1, pul_max
          ! diagonal elements of Aij
          RR = dot(n_my_grid_points, cov_resid_pul(:,ii,spin), 1, &
                   resid_pul(:,ii,spin), 1)
          call gsum(RR)
          Aij(ii,ii,spin) = RR
          ! the rest
          if (ii > 1) then
             do jj = 1, ii - 1
                RR = dot(n_my_grid_points, cov_resid_pul(:,ii,spin), 1, &
                         resid_pul(:,jj,spin), 1)
                call gsum(RR)
                Aij(ii,jj,spin) = RR
                ! Aij is symmetric even with wave-dependent metric
                Aij(jj,ii,spin) = RR
             end do ! jj
          end if ! (ii > 1)
       end do ! ii
    end do ! spin

    ! solve for alpha
    call DoPulay(iPulay, Aij, alpha, pul_max, maxpulaySC, inode, &
                 ionode)

    ! Compute the optimal rho and do mixing
    do spin = 1, nspin
       rho_opt(:,spin) = zero
       do ii = 1, pul_max
          call axpy(n_my_grid_points, alpha(ii,spin), &
                    rho_pul(:,ii,spin), 1, rho_opt(:,spin), 1)
          call axpy(n_my_grid_points, mix_factor(spin) * alpha(ii,spin), &
                    k_resid_pul(:,ii,spin), 1, rho_opt(:,spin), 1)
       end do
    end do

!****lat<$
    call stop_timer(t=tmr_std_loc,who='get_pulay_optimal_rho',echo=.true.)
!****lat>$

    return
  end subroutine get_pulay_optimal_rho
  !*****


  !!****f* SelfCon_module/update_pulay_history
  !! PURPOSE
  !!   From input rho_in, calculated the corresponding residual and
  !!   electron numbers and store in iPulay slot of the pulay history,
  !!   replacing the original content in the slot.
  !! USAGE
  !! INPUTS
  !! OUTPUT
  !! RETURN VALUE
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/05/22
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  subroutine update_pulay_history(iPulay, rho_in, reset_L,             &
                                  fixed_potential, vary_mu,            &
                                  n_L_iterations, L_tol, total_energy, &
                                  rho_pul, resid_pul, k_resid_pul,     &
                                  cov_resid_pul)
    use datatypes
    use GenBlas
    use dimens,         only: n_my_grid_points, grid_point_volume
    use EarlySCMod,     only: get_new_rho
    use GenComms,       only: gsum, cq_abort
    use hartree_module, only: kerker, kerker_and_wdmetric, wdmetric
    use global_module,  only: nspin, flag_fix_spin_population
    use maxima_module,  only: maxngrid
    use memory_module,  only: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    ! Passed parameters
    integer,      intent(in)  :: iPulay, n_L_iterations
    logical,      intent(in)  :: reset_L, fixed_potential, vary_mu
    real(double), intent(in)  :: L_tol
    real(double), intent(out) :: total_energy
    real(double), dimension(maxngrid,nspin),             intent(in)    :: &
         rho_in
    real(double), dimension(maxngrid,maxpulaySC,nspin),  intent(inout) :: &
         rho_pul, resid_pul, k_resid_pul, cov_resid_pul
    ! local variables
    integer :: spin, stat
    real(double), dimension(:,:), allocatable :: rho_out, resid
    type(cq_timer) :: tmr_std_loc

!****lat<$
    call start_timer(t=tmr_std_loc,who='update_pulay_history',where=5,&
         level=2,echo=.true.)
!****lat>$

    allocate(rho_out(maxngrid,nspin), resid(maxngrid,nspin), STAT=stat)
    if (stat /= 0) &
         call cq_abort("update_pulay_history: Error alloc mem: ", maxngrid, nspin)
    call reg_alloc_mem(area_SC, 2*maxngrid*nspin, type_dbl)

    ! store rho_in to history
    do spin = 1, nspin
       rho_pul(1:n_my_grid_points,iPulay,spin) = &
            rho_in(1:n_my_grid_points,spin)
    end do

    ! calculate residual of rho_in
    rho_out = zero
    resid = zero
    call get_new_rho(.false., reset_L, fixed_potential, vary_mu,   &
                     n_L_iterations, L_tol, total_energy, rho_in, &
                     rho_out, maxngrid)
    do spin = 1, nspin
       resid(1:n_my_grid_points,spin) = &
            rho_out(1:n_my_grid_points,spin) - &
            rho_in(1:n_my_grid_points,spin)
       ! replace the resid_pul history at iPulay to new value
       resid_pul(1:n_my_grid_points,iPulay,spin) = &
            resid(1:n_my_grid_points,spin)
    end do

    ! calculate the Kerker preconditioned or wave-dependent metric
    ! covarient residuals if requires to
    do spin = 1, nspin
       if (flag_Kerker) then
          if (flag_wdmetric) then
             call kerker_and_wdmetric(resid(:,spin), &
                                      cov_resid_pul(:,iPulay,spin), &
                                      maxngrid, q0, q1)
          else
             call kerker(resid(:,spin), maxngrid, q0)
          end if
          ! replace the k_resid_pul history at iPulay to new value
          k_resid_pul(1:n_my_grid_points,iPulay,spin) = &
               resid(1:n_my_grid_points,spin)
       else
          if (flag_wdmetric) then
             call wdmetric(resid(:,spin), cov_resid_pul(:,iPulay,spin), &
                           maxngrid, q1)
          end if
       end if ! flag_Kerker
    end do ! spin

    deallocate(rho_out, resid, STAT=stat)
    if (stat /= 0) &
         call cq_abort("update_pulay_history: Error dealloc mem")
    call reg_dealloc_mem(area_SC, 2*maxngrid*nspin, type_dbl)

!****lat<$
    call stop_timer(t=tmr_std_loc,who='update_pulay_history',echo=.true.)
!****lat>$

    return
  end subroutine update_pulay_history
  !*****


  !!****f* SelfCon_module/get_atomic_charge *
  !!
  !!  NAME
  !!   get_atomic_charge
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Computes and prints atomic charges (Mulliken analysis)
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   M. Todorovic
  !!  CREATION DATE
  !!   2008/04/02
  !!
  !!  MODIFICATION HISTORY
  !!   2011/05/25 17:00 dave/ast
  !!    Changed assignment of unit to Conquest standard
  !!   2011/12/07 L.Tong
  !!    - Added spin polarisation
  !!      note that only the total atomic charges are printed to file
  !!    - Added registrations for memory usage
  !!   2012/03/21 L.Tong
  !!   - Rewrote spin implementation
  subroutine get_atomic_charge()

    use datatypes
    use numbers
    use global_module,  only: ni_in_cell, area_SC, nspin, spin_factor
    use primary_module, only: bundle
    use matrix_data,    only: Srange
    use mult_module,    only: matK, matS, atom_trace, matrix_sum, &
                              allocate_temp_matrix, free_temp_matrix
    use atoms,          only: atoms_on_node
    use GenComms,       only: gsum, inode, ionode, cq_abort
    use input_module,   only: io_assign, io_close
    use memory_module,  only: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    ! Local variables
    integer :: chun, stat, n,l, glob_ind, spin
    integer, dimension(nspin) :: temp_mat
    real(double), dimension(:),   allocatable :: charge
    real(double), dimension(:,:), allocatable :: node_charge

    call start_timer(tmr_std_chargescf)
    ! prepare arrays and suitable matrices
    l = bundle%n_prim
    call start_timer(tmr_std_allocation)
    allocate(charge(ni_in_cell), STAT=stat)
    if (stat /= 0) call cq_abort("get_atomic_charge: Error alloc mem: ", ni_in_cell)
    call reg_alloc_mem(area_SC, ni_in_cell, type_dbl)
    allocate(node_charge(l,nspin), STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error allocating charge arrays in get_atomic_charge.", &
                       stat)
    call reg_alloc_mem(area_SC, nspin * l, type_dbl)
    call stop_timer(tmr_std_allocation)
    ! perform charge summation
    ! automatically called on each node
    node_charge = zero
    charge = zero
    do spin = 1, nspin
       temp_mat(spin) = allocate_temp_matrix(Srange, 0)
       call matrix_sum(zero, temp_mat(spin), one, matK(spin))
       call atom_trace(temp_mat(spin), matS, l, node_charge(:,spin))
       ! sum from the node_charge into the total charge array
       do n = 1, l
          glob_ind = atoms_on_node(n, inode)
          charge(glob_ind) = charge(glob_ind) + &
                             spin_factor * node_charge(n,spin)
       end do
    end do
    call gsum(charge, ni_in_cell)

    ! output
    if (inode == ionode) then
       ! write (*, *) 'Writing charge on individual atoms...'
       call io_assign(chun)
       open(unit = chun, file='AtomCharge.dat')
       do n = 1, ni_in_cell
          write (chun, fmt='(f15.10)') charge(n)
       end do
       call io_close(chun)
    end if

    call start_timer(tmr_std_allocation)
    deallocate(node_charge, STAT=stat)
    if (stat /= 0) &
         call cq_abort("Error deallocating charge arrays in get_atomic_charge.", &
                       stat)
    call reg_dealloc_mem(area_SC, nspin * l, type_dbl)
    deallocate(charge, STAT=stat)
    if (stat /= 0) call cq_abort("get_atomic_charge: Error dealloc mem")
    call reg_dealloc_mem(area_SC, ni_in_cell, type_dbl)
    do spin = nspin, 1, -1
       call free_temp_matrix(temp_mat(spin))
    end do
    call stop_timer (tmr_std_allocation)
    call stop_timer (tmr_std_chargescf)

    return
  end subroutine get_atomic_charge
  !!***


end module SelfCon

