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
!!   Contains master subroutines used in self-consistency, and 
!!   various useful variables.  Both "early" and "late" stage
!!   self-consistency strategies are described in the notes 
!! "How to achieve self-consistency" by MJG to be found in the
!!   Conquest documentation repository, also in Chem. Phys. Lett.
!!   325, 796 (2000) and in pieces in Cquest (my directory NewSC)
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
!!    Small changes: made residuals square roots of dot product (as they should be)
!!   10:35, 06/03/2003 drb 
!!    Added simple linear mixing and one or two other small changes
!!   2006/03/06 05:57 dave
!!    Rewrote calls to entire module to simplify
!!   2006/09/20 17:06 dave
!!    Made parameters user-definable
!!   2008/02/04 08:23 dave
!!    Changed for output to file not stdout
!!  SOURCE
!!
module SelfCon

  use datatypes
  use global_module, ONLY: io_lun

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

  real(double), save :: A
  real(double), save :: q0
  logical, save :: flag_linear_mixing
  real(double), save :: EndLinearMixing

  ! -------------------------------------------------------
  ! RCS ident string for object file id
  ! -------------------------------------------------------
  character(len=80), save, private :: RCSid = "$Id$"

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
!!  SOURCE
!!
  subroutine new_SC_potl( record, self_tol, reset_L, fixed_potential, vary_mu, n_L_iterations, &
       number_of_bands, L_tol, mu, total_energy)

    use datatypes
    use PosTan, ONLY: PulayC, PulayBeta, SCC, SCBeta, pos_tan, &
         max_iters, SCE, SCR, fit_coeff
    use numbers
    use global_module, ONLY: iprint_SC, flag_self_consistent, flag_SCconverged
    use H_matrix_module, ONLY: get_H_matrix
    use DMMin, ONLY: FindMinDM
    use energy, ONLY: get_energy
    use GenComms, ONLY: inode, ionode, cq_abort
    use dimens, ONLY: n_my_grid_points
    use density_module, ONLY: density
    use maxima_module, ONLY: maxngrid

    implicit none

    ! Passed variables
    logical :: record   ! Flags whether to record dE vs R
    logical vary_mu, fixed_potential, reset_L

    integer n_L_iterations

    real(double) :: self_tol
    real(double) ::  number_of_bands, L_tol, mu
    real(double) :: total_energy

    ! Local variables
    integer :: mod_early, ndone, i, nkeep, ndelta, stat
    logical :: done, problem, early
    real(double) :: SC_tol, DMM_tol, LastE, electrons

    ! Build H matrix *with NL and KE*
    call get_H_matrix(.true., fixed_potential, electrons, density, maxngrid)
    if(.NOT.flag_self_consistent) then
       call FindMinDM(n_L_iterations, number_of_bands, vary_mu, L_tol, mu, inode, ionode, reset_L, .false.)
       call get_energy(total_energy)
       return
    end if
    if(inode==ionode) write(io_lun,fmt='(8x,"Starting self-consistency.  Tolerance: ",e12.5,/)') self_tol
    if(record) then
       if(inode==ionode.AND.iprint_SC>1) write(io_lun,*) 'Original tol: ',L_tol
       SC_tol = self_tol
       DMM_tol = L_tol
    else
       SC_tol = max(self_tol,SCC*(self_tol**SCBeta))
       DMM_tol = max(L_tol,PulayC*((0.1_double*self_tol)**PulayBeta))
    endif
    ndone = 0
    done = .false.
    problem = .false.
    early = .false.
    flag_SCconverged = .true.
    ! Check on whether we need to do early iterations
    if(.NOT.allocated(EarlyRecord)) then
       allocate(EarlyRecord(maxearlySC),STAT=stat)
       if(stat/=0) call cq_abort("Error allocating EarlyRecord in SCpotl: ",maxearlySC)
       EarlyRecord = 0
    end if
    do i=1,maxearlySC
       if(EarlyRecord(i)==0) early = .true.
       if(inode.eq.ionode.AND.iprint_SC>1) write(io_lun,*) 'early: ',i,EarlyRecord(i)
    enddo
    ! Loop until self-consistent
    do while(.NOT.done)
       if(flag_linear_mixing) then
          !call LinearMixSC(done,ndone, SC_tol, &
          call PulayMixSCA(done,ndone, SC_tol, reset_L, fixed_potential, vary_mu, n_L_iterations, &
               number_of_bands, DMM_tol, mu, total_energy, density, maxngrid)
       else if(early.OR.problem) then ! Early stage strategy
          earlyL = .false.
          reset_L = .true.
          if(inode==ionode.AND.iprint_SC>0) write(io_lun,*) '********** EarlySC **********'
          if(inode==ionode.AND.problem) write(io_lun,*) 'Problem restart'
          call earlySC(record,done,earlyL,ndone, SC_tol, reset_L, fixed_potential, vary_mu, n_L_iterations, &
               number_of_bands, DMM_tol, mu, total_energy, density, maxngrid)
          ! Check for early/late stage switching and update variables
          mod_early = 1+mod(earlyIters,maxearlySC)
          if(earlyL) EarlyRecord(mod_early) = 1
       endif
       if(.NOT.done) then ! Late stage strategy
          reset_L = .true. !.false.
          if(inode==ionode.AND.iprint_SC>0) write(io_lun,*) '********** LateSC **********'
          call lateSC(record,done,ndone, SC_tol, reset_L, fixed_potential, vary_mu, n_L_iterations, &
               number_of_bands, DMM_tol, mu, total_energy, density, maxngrid)
          problem = .true.  ! This catches returns from lateSC when not done
       endif
       ! Now decide whether we need to bother going to early
       early = .false.
       do i=1,maxearlySC
          if(EarlyRecord(i)==0) early = .true.
          if(inode.eq.ionode.AND.iprint_SC>1) write(io_lun,*) 'early: ',i,EarlyRecord(i)
       enddo
    enddo
    if(record) then ! Fit the C and beta coefficients
       if(inode==ionode.AND.iprint_SC>1) then
          write(io_lun,*) '  List of residuals and energies'
          if(ndone>0) then
             do i=1,ndone
                write(io_lun,7) i, SCR(i), SCE(i)
             enddo
          end if
       endif
       call fit_coeff(SCC, SCBeta, SCE, SCR, ndone)
       if(inode==ionode.AND.iprint_SC>1) write(io_lun,6) SCC,SCBeta
    endif
6   format(8x,'dE to dR parameters - C: ',f15.8,' beta: ',f15.8)
7   format(8x,i4,2f15.8)
    return
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
!!  SOURCE
!!
  subroutine earlySC(record,done,EarlyLin,ndone,self_tol, reset_L, fixed_potential, vary_mu, n_L_iterations, &
       number_of_bands, L_tol, mu, total_energy, rho, size)

    use datatypes
    use numbers
    use GenBlas
    use PosTan
    use EarlySCMod
    use GenComms, ONLY: cq_abort, gsum, my_barrier, inode, ionode
    use dimens, ONLY: n_my_grid_points, grid_point_volume
    use global_module, ONLY: ne_in_cell, area_SC, flag_continue_on_SC_fail, flag_SCconverged
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use maxima_module, ONLY: maxngrid

    implicit none

    ! Passed variables
    logical :: record, done, EarlyLin
    logical :: vary_mu, fixed_potential, reset_L

    integer :: ndone, size
    integer :: n_L_iterations

    real(double), dimension(size) :: rho
    real(double) :: self_tol
    real(double) :: number_of_bands, L_tol, mixing, mu
    real(double) :: total_energy

    ! Local variables
    integer :: n_iters,irc,ierr
    logical :: linear, left, moved, init_reset, Lrec
    real(double) :: R0,R1,Rcross, Rcross_a, Rcross_b, Rcross_c
    real(double) :: Ra, Rb, Rc, E0, E1, dE, stLtol
    real(double) :: lambda_1, lambda_a, lambda_b, lambda_c
    real(double) :: zeta_exact, ratio, qatio, lambda_up
    real(double) :: pred_Rb, Rup, Rcrossup, zeta_num, zeta
    real(double), allocatable, dimension(:) :: resid0, rho1, residb

    allocate(resid0(maxngrid),rho1(maxngrid),residb(maxngrid),STAT=ierr)
    if(ierr/=0) call cq_abort("Allocation error in earlySC: ",maxngrid, ierr)
    call reg_alloc_mem(area_SC, 3*maxngrid, type_dbl)
    !    write(io_lun,*) inode,'Welcome to earlySC'
    init_reset = reset_L
    call my_barrier()
    ! Test to check for available iterations
    n_iters = ndone
    if(n_iters>=maxitersSC) then
       if(.NOT.flag_continue_on_SC_fail) call cq_abort('earlySC: Too many self-consisteny iterations: ',n_iters,&
            maxitersSC)
       flag_SCconverged = .false.
       done = .true.
       return
    endif
    ! Decide on whether or not to record dE vs R for L min
    if(record) then 
       Lrec = .true.
    else
       Lrec = .false.
    endif
    ! Compute residual of initial density
    call get_new_rho(Lrec, reset_L, fixed_potential, vary_mu, n_L_iterations, number_of_bands, L_tol, mu, &
         total_energy, rho, rho1, maxngrid)
    Lrec = .false.
    resid0 = rho1 - rho
    R0 = dot(n_my_grid_points,resid0,1,resid0,1)
    call gsum(R0)
    R0 = sqrt(grid_point_volume*R0)/ne_in_cell
    if(inode.eq.ionode) write(io_lun,*) 'In EarlySC, R0 is ',R0
    if(R0<self_tol) then
       done = .true.
       if(inode==ionode) write(io_lun,*) 'Done ! Self-consistent'
    endif
    linear = .false.
    lambda_1 = InitialLambda
    ! Now loop until we achieve a reasonable reduction or linearity
    ! -------------------------------------------------------------
    ! START OF MAIN LOOP
    ! -------------------------------------------------------------
    E0 = total_energy
    dE = L_tol
    do while((.NOT.done).AND.(.NOT.linear).AND.(n_iters<maxitersSC))
       if(init_reset) reset_L = .true.
       n_iters = n_iters + 1
       ! Go to initial lambda and find R
       if(inode==ionode) write(io_lun,*) '********** Early iter ',n_iters
       R1 = getR2( MixLin, lambda_1, reset_L, fixed_potential, vary_mu, n_L_iterations, number_of_bands, L_tol, mu, &
            total_energy, rho, rho1,residb, maxngrid)
       Rcross = dot(n_my_grid_points, residb, 1, resid0, 1)
       if(R0<self_tol) then
          done = .true.
          if(inode==ionode) write(io_lun,*) 'Done ! Self-consistent'
          exit
       endif
       ! Test for acceptable reduction
       if(inode.eq.ionode) write(io_lun,*) 'In EarlySC, R1 is ',R1
       if(R1<=thresh*R0) then
          if (inode==ionode) write(io_lun,*) 'No need to reduce'
       else
          if (inode==ionode) write(io_lun,*) 'Reducing lambda'
          call reduceLambda(R0,R1,lambda_1, thresh, MixLin, reset_L, fixed_potential, vary_mu, n_L_iterations, &
               maxitersSC, number_of_bands, L_tol, mu, total_energy,rho,rho1,residb, maxngrid)
          Rcross = dot(n_my_grid_points, residb, 1, resid0, 1)
          call gsum(Rcross)
       endif ! Preliminary reduction of lambda
       ! At this point, we have two points: 0 and R0, lambda_1 and R1
       ! Now search for a good bracket
       if(R1>R0) then  ! Take (a,b) = (lam,0)
          left = .false.
          lambda_a = lambda_1
          Ra = R1
          Rcross_a = Rcross
          lambda_b = zero
          Rb = R0
          Rcross_b = R0
          call copy(n_my_grid_points, resid0, 1, residb, 1)
       else            ! Take (a,b) = (0,lam)
          left = .true.
          lambda_a = zero
          Ra = R0
          Rcross_a = R0
          lambda_b = lambda_1
          Rb = R1
          Rcross_b = Rcross
          moved = .true.
       endif
       call bracketMin(moved, &
            lambda_a, lambda_b, lambda_c, Rcross_a, Rcross_b, Rcross_c, &
            Ra, Rb, Rc, MixLin, reset_L, fixed_potential, vary_mu, &
            n_L_iterations, maxitersSC, number_of_bands, L_tol, mu,&
            total_energy, rho, rho1, residb, resid0, maxngrid)
       ! Test for acceptability - otherwise do line search
       reset_L = .false.
       if(inode==ionode) write(io_lun,*) 'a,b,c: ',lambda_a, lambda_b, lambda_c, &
            Ra,Rb,Rc
       if(moved.AND.Rb<ReduceLimit*R0) then ! Accept
          if (inode==ionode) &
               write(io_lun,*) 'Line search accepted after bracketing', R0, Rb
       else ! Line search
          if(.NOT.moved) then
             if (inode==ionode) &
                  write(io_lun,*) 'Search because intern point not moved'
          else
             if (inode==ionode) &
                  write(io_lun,*) 'Search because R reduction too small'
          endif
          ! Do golden section or brent search
          ! call searchMin(lambda_a, lambda_b, lambda_c, Rcross_a, Rcross_c, &
          call brentMin(lambda_a, lambda_b, lambda_c, Rcross_a, Rcross_c, &
               Ra, Rb, Rc, MixLin, reset_L, fixed_potential, vary_mu, &
               n_L_iterations, maxitersSC, number_of_bands, L_tol, mu,&
               total_energy, rho, rho1, residb, maxngrid)
       endif ! end if(reduced)
       if(abs(lambda_a)>abs(lambda_c)) then 
          Rcrossup = Rcross_a
          lambda_up = lambda_a
          Rup = Ra
       else
          Rcrossup = Rcross_c
          lambda_up = lambda_c
          Rup = Rc
       endif
       E1 = total_energy
       dE = E1 - E0
       E0 = E1
       ! Now test for linearity
       zeta_exact = abs(Rb - R0)
       ratio = lambda_b/lambda_up
       qatio = one - ratio
       pred_Rb = qatio*qatio*R0+ratio*ratio*Rup+two*ratio*qatio*Rcrossup
       zeta_num = abs(Rb - pred_Rb)
       zeta = zeta_num/zeta_exact
       if (inode==ionode) write(io_lun,31) zeta
31     format('Linearity monitor for this line search: ',e15.6)
       if(zeta<crit_lin) then
          linear = .true.
          if(inode==ionode) write(io_lun,*) 'Linearity fulfilled'
       endif
       ! And record residual and total energy if necessary
       if(record) then
          SCE(n_iters) = total_energy
          SCR(n_iters) = Rb
       endif
       resid0(1:n_my_grid_points) = residb(1:n_my_grid_points)
       call mixtwo(n_my_grid_points, MixLin, lambda_b, &
            rho, rho1, residb)
       rho = zero
       rho(1:n_my_grid_points) = residb(1:n_my_grid_points)
       rho1 = zero
       rho1(1:n_my_grid_points) = rho(1:n_my_grid_points) + &
            resid0(1:n_my_grid_points)
       lambda_1 = lambda_b
       R0 = dot(n_my_grid_points,resid0,1,resid0,1)
       call gsum(R0)
!       R0 = sqrt(R0)
       if(inode.eq.ionode) write(io_lun,*) 'In EarlySC, R0 is ',R0
       if(R0<self_tol) then
          done = .true.
          if(inode==ionode) write(io_lun,*) 'Done ! Self-consistent'
       endif
    enddo ! end do while(.NOT.linear.AND..NOT.done)
    ! -------------------------------------------------------------
    ! END OF MAIN LOOP
    ! -------------------------------------------------------------
    ! This helps us to decide when to switch from early to late automatically
    if(linear.AND.n_iters==1) EarlyLin = .true.
    ndone = n_iters
    deallocate(resid0,rho1,residb,STAT=ierr)
    if(ierr/=0) call cq_abort("Deallocation error in earlySC: ",maxngrid, ierr)
    call reg_dealloc_mem(area_SC, 3*maxngrid, type_dbl)
    return
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
!!  SOURCE
!!
  subroutine lateSC(record,done,ndone,self_tol, reset_L, fixed_potential, vary_mu, n_L_iterations, &
       number_of_bands, L_tol, mu, total_energy,rho, size)

    use datatypes
    use numbers
    use PosTan
    use Pulay, ONLY: DoPulay2D
    use GenBlas, ONLY: dot, rsum
    use dimens, ONLY: n_my_grid_points, grid_point_volume
    use EarlySCMod, ONLY: get_new_rho
    use GenComms, ONLY: gsum, cq_abort, inode, ionode
    use hartree_module, ONLY: kerker
    use global_module, ONLY: ne_in_cell, area_SC, flag_continue_on_SC_fail, flag_SCconverged
    use maxima_module, ONLY: maxngrid
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    ! Passed variables
    logical :: record, done
    logical :: vary_mu, fixed_potential, reset_L

    integer :: size
    integer :: ndone
    integer :: n_L_iterations

    real(double) :: self_tol
    real(double) :: number_of_bands, L_tol, mixing, mu
    real(double) :: total_energy
    real(double), dimension(size) :: rho

    ! Local variables
    integer :: n_iters, n_pulay, npmod, pul_mx, i, j, stat
    logical :: linear
    real(double), allocatable, dimension(:,:) :: rho_pul
    real(double), allocatable, dimension(:,:) :: resid_pul
    real(double), allocatable, dimension(:) :: rho1
    real(double), allocatable, dimension(:) :: Kresid
    real(double), allocatable, dimension(:,:) :: Aij
    real(double), allocatable, dimension(:) :: alph
    real(double) :: R0, R, R1, E0, E1, dE, stLtol, tmp

    allocate(rho_pul(maxngrid,maxpulaySC), resid_pul(maxngrid,maxpulaySC),rho1(maxngrid), &
         Kresid(maxngrid),Aij(maxpulaySC, maxpulaySC), alph(maxpulaySC), STAT=stat)
    if(stat/=0) call cq_abort("Allocation error in lateSC: ",maxngrid, maxpulaySC)
    call reg_alloc_mem(area_SC, 2*maxngrid*(1+maxpulaySC)+maxpulaySC*(maxpulaySC+1), type_dbl)
    done = .false.
    linear = .true.
    n_iters = ndone
    if(n_iters>=maxitersSC) then
       if(.NOT.flag_continue_on_SC_fail) call cq_abort('lateSC: too many iterations: ',n_iters, maxitersSC)
       flag_SCconverged = .false.
       done = .true.
       return
    endif
    ! Compute residual of initial density
    rho_pul(1:n_my_grid_points, 1) = rho(1:n_my_grid_points)
    call get_new_rho(.false., reset_L, fixed_potential, vary_mu, n_L_iterations, number_of_bands, L_tol, mu, &
         total_energy, rho, rho1, maxngrid)
    resid_pul(1:n_my_grid_points,1) = rho1(1:n_my_grid_points) - &
         rho(1:n_my_grid_points)
    !Kresid(1:n_my_grid_points) = resid_pul(1:n_my_grid_points,1)
    !call kerker(Kresid,maxngrid,q0)
    call kerker(resid_pul(:,1),maxngrid,q0)
    R0 = dot(n_my_grid_points,resid_pul(:,1),1,resid_pul(:,1),1)
    call gsum(R0)
    R0 = sqrt(grid_point_volume*R0)/ne_in_cell
    ! Try this ?!
    if(inode==ionode) write(io_lun,*) 'Residual is ',R0
    ! Old output becomes new input
    rho(1:n_my_grid_points) = rho(1:n_my_grid_points) + resid_pul(1:n_my_grid_points,1)
    tmp = grid_point_volume*rsum(n_my_grid_points,rho,1)
    call gsum(tmp)
    rho = ne_in_cell*rho/tmp
    ! Now loop until we achieve SC or something goes wrong
    n_pulay = 0
    do while((.NOT.done).AND.(linear).AND.(n_iters<maxitersSC))
       if(R0<1.0_double) reset_L = .false.
       n_iters = n_iters + 1
       if(inode==ionode) write(io_lun,*) '********** Late iter ',n_iters
       n_pulay = n_pulay+1
       ! Storage for pulay charges/residuals
       npmod = mod(n_pulay, maxpulaySC)+1
       pul_mx = min(n_pulay+1, maxpulaySC)
       if(inode==ionode) write(io_lun,*) 'npmod, pul_mx: ',npmod,pul_mx
       ! For the present output, find the residual (i.e. R_n^\prime)
       call get_new_rho(.false., reset_L, fixed_potential, vary_mu, n_L_iterations, number_of_bands, L_tol, mu, &
            total_energy, rho, rho1, maxngrid)
       rho_pul(1:n_my_grid_points, npmod) = rho(1:n_my_grid_points)
       resid_pul(1:n_my_grid_points, npmod) = rho1(1:n_my_grid_points) - &
            rho(1:n_my_grid_points)
       !Kresid(1:n_my_grid_points) = resid_pul(1:n_my_grid_points,npmod)
       !call kerker(Kresid,maxngrid,q0)
       call kerker(resid_pul(:,npmod),maxngrid,q0)
       ! Do the pulay-style mixing
       ! Form the A matrix (=<R_i|R_j>)
       do i=1,pul_mx
          do j=1,pul_mx
             R = dot(n_my_grid_points, resid_pul(:,i), 1, resid_pul(:,j),1)
             call gsum(R)
             Aij(j,i) = R
          enddo
       enddo
       !      if(inode==ionode)write(io_lun,*) 'Aij: ',Aij(1:pul_mx,1:pul_mx)
       ! Solve to get alphas
       call DoPulay2D(Aij,alph,pul_mx,maxpulaySC,inode,ionode)
       if(inode==ionode)write(io_lun,*) 'alph: ',alph(1:pul_mx)
       ! Build new input density - we could do this wrap around, but will 
       ! need rho anyway, so might as well use it.
       rho = zero
       do i=1,pul_mx
          rho(1:n_my_grid_points) = rho(1:n_my_grid_points) +  &
               alph(i)*rho_pul(1:n_my_grid_points,i)
       enddo
       tmp = grid_point_volume*rsum(n_my_grid_points,rho,1)
       call gsum(tmp)
       rho = ne_in_cell*rho/tmp
       rho_pul(1:n_my_grid_points,npmod) = rho(1:n_my_grid_points)
       ! Generate the residual either exactly or by extrapolation
       if(mod(n_iters, n_exact)==0) then
          if(inode==ionode) write(io_lun,*) 'Generating rho exactly'
          resid_pul(1:n_my_grid_points, npmod) = alph(npmod)* &
               resid_pul(1:n_my_grid_points, npmod)
          if(pul_mx==maxpulaySC) then
             do i=npmod+1, npmod+maxpulaySC-1
                j = mod(i,maxpulaySC)
                if(j==0) j=maxpulaySC
                resid_pul(1:n_my_grid_points,npmod) = &
                     resid_pul(1:n_my_grid_points, npmod) + &
                     alph(j)*resid_pul(1:n_my_grid_points, j)
             enddo
          else
             do i=1,npmod-1 
                resid_pul(1:n_my_grid_points,npmod) = &
                     resid_pul(1:n_my_grid_points, npmod) + &
                     alph(i)*resid_pul(1:n_my_grid_points, i)
             enddo
          endif
          R1 = dot(n_my_grid_points, resid_pul(:,npmod),1,resid_pul(:,npmod),1)
          call gsum(R1)
          R1 = sqrt(grid_point_volume*R1)/ne_in_cell
          if(inode==ionode) write(io_lun,fmt='(8x,"Predicted residual is ",f20.12)') R1
          call get_new_rho(.false., reset_L, fixed_potential, vary_mu, n_L_iterations, number_of_bands, L_tol, mu, &
               total_energy, rho, rho1, maxngrid)
          resid_pul(1:n_my_grid_points, npmod) = rho1(1:n_my_grid_points) - &
               rho(1:n_my_grid_points)
          !Kresid(1:n_my_grid_points) = resid_pul(1:n_my_grid_points,npmod)
          !call kerker(Kresid,maxngrid,q0)
          call kerker(resid_pul(:,npmod),maxngrid,q0)
          R1 = dot(n_my_grid_points, resid_pul(:,npmod),1,resid_pul(:,npmod),1)
          call gsum(R1)
          R1 = sqrt(grid_point_volume*R1)/ne_in_cell
          if(inode==ionode) write(io_lun,fmt='(8x,"Actual residual is ",f20.12)') R1
          rho(1:n_my_grid_points) = rho(1:n_my_grid_points) + resid_pul(1:n_my_grid_points,npmod)
          tmp = grid_point_volume*rsum(n_my_grid_points,rho,1)
          call gsum(tmp)
          rho = ne_in_cell*rho/tmp
       else
          if(inode==ionode) write(io_lun,*) 'Generating rho by interpolation'
          ! Clever, wrap-around way to do it
          resid_pul(1:n_my_grid_points, npmod) = alph(npmod)* &
               resid_pul(1:n_my_grid_points, npmod)
          if(pul_mx==maxpulaySC) then
             do i=npmod+1, npmod+maxpulaySC-1
                j = mod(i,maxpulaySC)
                if(j==0) j=maxpulaySC
                resid_pul(1:n_my_grid_points,npmod) = &
                     resid_pul(1:n_my_grid_points, npmod) + &
                     alph(j)*resid_pul(1:n_my_grid_points, j)
             enddo
          else
             do i=1,npmod-1 
                resid_pul(1:n_my_grid_points,npmod) = &
                     resid_pul(1:n_my_grid_points, npmod) + &
                     alph(i)*resid_pul(1:n_my_grid_points, i)
             enddo
          endif
          ! Output for free if we've extrapolated
          rho(1:n_my_grid_points) = rho_pul(1:n_my_grid_points, npmod) + &
               resid_pul(1:n_my_grid_points, npmod)
          tmp = grid_point_volume*rsum(n_my_grid_points,rho,1)
          call gsum(tmp)
          rho = ne_in_cell*rho/tmp
       endif
       ! Test for (a) convergence and (b) linearity
       E1 = total_energy
       dE = E1 - E0
       E0 = E1
       R1 = dot(n_my_grid_points, resid_pul(:,npmod),1,resid_pul(:,npmod),1)
       call gsum(R1)
       R1 = sqrt(grid_point_volume*R1)/ne_in_cell
       if(R1>2.0_double*R0.AND.mod(n_iters, n_exact)/=0) then
          linear = .false.
          if(npmod>1) then
             rho(1:n_my_grid_points) = rho_pul(1:n_my_grid_points, npmod-1)
          else
             rho(1:n_my_grid_points) = rho_pul(1:n_my_grid_points, pul_mx)
          endif
          write(io_lun,*) 'PANIC ! Residual increase !'
       endif
       R0 = R1
       if(R0<self_tol) then
          done = .true.
          if(inode==ionode) write(io_lun,*) 'Done ! Self-consistent'
       endif
       if(inode==ionode) write(io_lun,*) 'Residual is ',R0
       ! Output an energy
       !      write(io_lun,*) 'Energy is: '
       if(record) then
          SCE(n_iters) = total_energy
          SCR(n_iters) = R0
       endif
    enddo
    ndone = n_iters
    if(inode==ionode) write(io_lun,*) 'Finishing lateSC after ',ndone, &
         'iterations with residual of ',R0
    deallocate(rho_pul, Kresid, resid_pul,rho1, Aij, alph, STAT=stat)
    if(stat/=0) call cq_abort("Dellocation error in lateSC: ",maxngrid, maxpulaySC)
    call reg_dealloc_mem(area_SC, 2*maxngrid*(1+maxpulaySC)+maxpulaySC*(maxpulaySC+1), type_dbl)
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
!!  SOURCE
!!
  subroutine LinearMixSC(done,ndone,self_tol, reset_L, fixed_potential, vary_mu, n_L_iterations, &
       number_of_bands, L_tol, mu, total_energy, rho, size)
    use datatypes
    use numbers
    use PosTan
    use Pulay, ONLY: DoPulay2D
    use GenBlas
    use dimens, ONLY: n_my_grid_points, grid_point_volume
    use EarlySCMod, ONLY: get_new_rho
    use GenComms, ONLY: gsum, cq_abort, inode, ionode
    use io_module, ONLY: dump_charge
    use hartree_module, ONLY: kerker
    use global_module, ONLY: ne_in_cell, area_SC, flag_continue_on_SC_fail, flag_SCconverged
    use maxima_module, ONLY: maxngrid
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    ! Passed variables
    logical :: done
    logical :: vary_mu, fixed_potential, reset_L

    integer :: ndone,size
    integer :: n_L_iterations

    real(double) :: self_tol
    real(double) :: number_of_bands, L_tol, mixing, mu
    real(double) :: total_energy
    real(double), dimension(size) :: rho
    ! Local variables
    integer :: n_iters, n_pulay, npmod, pul_mx, i, j, stat
    real(double), allocatable, dimension(:) :: rho1, resid
    real(double) :: R0, R, R1, E0, E1, dE, stLtol

    allocate(rho1(maxngrid), resid(maxngrid), STAT=stat)
    if(stat/=0) call cq_abort("Error allocating rho1 in linearMixSC: ",maxngrid,stat)
    call reg_alloc_mem(area_SC, 2*maxngrid, type_dbl)
    done = .false.
    n_iters = ndone
    if(n_iters>=maxitersSC) then
       if(.NOT.flag_continue_on_SC_fail) call cq_abort('lateSC: too many iterations: ',n_iters, maxitersSC)
       flag_SCconverged = .false.
       done = .true.
       return
    endif
    R0 = 100.0_double
    do while(R0>EndLinearMixing)
       ! Generate new charge and find residual
       call get_new_rho(.false., reset_L, fixed_potential, vary_mu, n_L_iterations, number_of_bands, L_tol, mu, &
            total_energy, rho, rho1, maxngrid)
       resid = rho1 - rho
       call kerker(resid,maxngrid,q0)
       R0 = dot(n_my_grid_points,resid,1,resid,1)
       call gsum(R0)
       R0 = sqrt(grid_point_volume*R0)/ne_in_cell
       if(inode==ionode) write(io_lun,*) 'Residual is ',R0
       ! Old output becomes new input
       !rho = A*rho1 + (1.0_double-A)*rho
       rho = rho + A*resid
       n_iters = n_iters+1
       call dump_charge(rho,size,inode)
    end do
    deallocate(rho1, resid, STAT=stat)
    if(stat/=0) call cq_abort("Error allocating rho1 in linearMixSC: ",maxngrid,stat)
    call reg_dealloc_mem(area_SC, 2*maxngrid, type_dbl)
!    ndone = n_iters
  end subroutine LinearMixSC
!!***

! VASP version
  subroutine PulayMixSC(done,ndone,self_tol, reset_L, fixed_potential, vary_mu, n_L_iterations, &
       number_of_bands, L_tol, mu, total_energy,rho,size)
    use datatypes
    use numbers
    use PosTan
    use Pulay, ONLY: DoPulay2D
    use GenBlas
    use dimens, ONLY: n_my_grid_points, grid_point_volume
    use EarlySCMod, ONLY: get_new_rho
    use GenComms, ONLY: gsum, cq_abort, inode, ionode
    use io_module, ONLY: dump_charge
    use hartree_module, ONLY: kerker
    use global_module, ONLY: ne_in_cell, area_SC,flag_continue_on_SC_fail, flag_SCconverged
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use maxima_module, ONLY: maxngrid

    implicit none

    ! Passed variables
    logical :: done
    logical :: vary_mu, fixed_potential, reset_L

    integer :: size
    integer :: ndone
    integer :: n_L_iterations

    real(double) :: self_tol
    real(double) :: number_of_bands, L_tol, mixing, mu
    real(double) :: total_energy
    real(double), dimension(size) :: rho

    ! Local variables
    integer :: n_iters, n_pulay, npmod, pul_mx, i, j, ii, m, next_i, stat
    real(double) :: R0, R, R1, E0, E1, dE, stLtol, norm
    real(double), allocatable, dimension(:,:) :: rho_pul, delta_rho, delta_R
    real(double), allocatable, dimension(:,:) :: resid_pul
    real(double), allocatable, dimension(:) :: rho1, resid, Kresid, old_resid
    real(double), allocatable, dimension(:,:) :: Aij
    real(double), allocatable, dimension(:) :: alph

    allocate(rho_pul(maxngrid,maxpulaySC), resid_pul(maxngrid,maxpulaySC),rho1(maxngrid), &
         delta_rho(maxngrid,maxpulaySC), delta_R(maxngrid,maxpulaySC), &
         resid(maxngrid), Kresid(maxngrid),old_resid(maxngrid), &
         Aij(maxpulaySC, maxpulaySC), alph(maxpulaySC), STAT=stat)
    if(stat/=0) call cq_abort("Allocation error in PulayMixSC: ",maxngrid, maxpulaySC)
    call reg_alloc_mem(area_SC, 4*maxngrid*(1+maxpulaySC)+maxpulaySC*(maxpulaySC+1), type_dbl)
    delta_R = zero
    delta_rho = zero
    done = .false.
    n_iters = ndone
    if(n_iters>=maxitersSC) then
       if(.NOT.flag_continue_on_SC_fail) call cq_abort('lateSC: too many iterations: ',n_iters, maxitersSC)
       flag_SCconverged = .false.
       done = .true.
       return
    endif
    ! m=1
    call get_new_rho(.false., reset_L, fixed_potential, vary_mu, n_L_iterations, number_of_bands, L_tol, mu, &
         total_energy, rho, rho1, size)
    ! Evaluate residual
    resid = rho1 - rho
    !Kresid = resid
    !call kerker(Kresid,maxngrid,q0)
    call kerker(resid,maxngrid,q0)
    R0 = dot(n_my_grid_points,resid,1,resid,1)
    call gsum(R0)
    R0 = sqrt(grid_point_volume*R0)/ne_in_cell
    if(inode==ionode) write(io_lun,*) 'Residual is ',R0
    ! Create new input charge
    !rho = rho + A*Kresid
    rho = rho + A*resid
    ! Store deltarho
    delta_rho = zero
    call axpy(n_my_grid_points,A,resid,1,delta_rho(:,1),1)
    !call axpy(n_my_grid_points,A,Kresid,1,delta_rho(:,1),1)
    n_iters = n_iters+1
    call dump_charge(rho,size,inode)    
    do m=2,maxitersSC
       ! calculate i (cyclical index for storing history) and pul_mx
       i = mod(m-2, maxpulaySC)+1
       next_i = mod(i,maxpulaySC)+1
       pul_mx = min(m-1, maxpulaySC)
       if(inode==ionode) write(io_lun,fmt='(8x,"Pulay iteration ",i4," counters are ",2i4)') m, i, pul_mx
       ! Generate new charge and find residual
       call get_new_rho(.false., reset_L, fixed_potential, vary_mu, n_L_iterations, number_of_bands, L_tol, mu, &
            total_energy, rho, rho1, size)
       call dump_charge(rho,size,inode)
       old_resid = resid
       resid = rho1 - rho ! F(m)
       call kerker(resid,n_my_grid_points,q0)
       R0 = dot(n_my_grid_points,resid,1,resid,1)
       call gsum(R0)
       R0 = sqrt(grid_point_volume*R0)/ne_in_cell
       if(inode==ionode) write(io_lun,*) 'Residual is ',R0
       if(R0<self_tol) then
          if(inode==ionode) write(io_lun,'(8x,"Reached tolerance")')
          done = .true.
          call dealloc_PulayMixSC
          return
       end if
       if(R0<EndLinearMixing) then
          if(inode==ionode) write(io_lun,'(8x,"Reached transition to LateSC")')
          call dealloc_PulayMixSC
          return
       end if
       ! Build deltaR(i)
       delta_R(:,i) = zero
       call axpy(n_my_grid_points,one,resid,1,delta_R(:,i),1) 
       call axpy(n_my_grid_points,-one,old_resid,1,delta_R(:,i),1)
       !call kerker(resid,n_my_grid_points,q0)
       ! Normalise
       !norm = dot(n_my_grid_points,delta_R(:,i),1,delta_R(:,i),1)
       !call gsum(norm)
       !norm = sqrt(norm)
       !delta_R(:,i) = delta_R(:,i)/norm
       !delta_rho(:,i) = delta_rho(:,i)/norm
       !if(inode==ionode) write(io_lun,fmt='(8x,"Norm of deltaF is ",f20.12)') norm
       ! now build new rho
       do ii=1,pul_mx
          do j=1,pul_mx
             R = dot(n_my_grid_points, delta_R(:,ii), 1, delta_R(:,j),1)
             call gsum(R)
             Aij(j,ii) = R
             if(ii==j) Aij(j,ii) = Aij(j,ii) !+ 0.0001_double ! w_0 = 0.01
          enddo
       enddo
       write(io_lun,*) 'A is ',Aij
       ! Solve to get alphas
       call DoPulay2D(Aij,alph,pul_mx,maxpulaySC,inode,ionode)
       write(io_lun,*) 'A is ',Aij
       alph = zero
       do j=1,pul_mx
          R = dot(n_my_grid_points, delta_R(:,j),1,resid(:),1)
          call gsum(R)
          do ii=1,pul_mx
             alph(ii) = alph(ii) + Aij(j,ii)*R
          enddo
       enddo
       if(inode==ionode)write(io_lun,*) 'alph: ',alph(1:pul_mx)
       ! Create delta rho - add on new rho
       !Kresid = resid
       !call kerker(Kresid,maxngrid,q0)
       !rho = rho + A*Kresid
       !delta_rho(:,next_i) = A*Kresid
       rho = rho + A*resid
       delta_rho(:,next_i) = A*resid
       do ii=1,pul_mx
          !Kresid = delta_R(:,ii)
          !call kerker(Kresid,maxngrid,q0)
          !rho = rho - alph(ii)*(delta_rho(:,ii)+A*Kresid)
          !delta_rho(:,next_i) = delta_rho(:,next_i) - alph(ii)*(delta_rho(:,ii)+A*Kresid)
          rho = rho - alph(ii)*(delta_rho(:,ii)+A*delta_R(:,ii))
          delta_rho(:,next_i) = delta_rho(:,next_i) - alph(ii)*(delta_rho(:,ii)+A*delta_R(:,ii))
       end do
       n_iters = n_iters+1
       call dump_charge(rho,size,inode)
    end do
    call dealloc_PulayMixSC
!    ndone = n_iters
    return

  contains

    subroutine dealloc_PulayMixSC
      
      deallocate(rho_pul, resid_pul,rho1, delta_rho, delta_R, resid, Kresid,old_resid, Aij, alph, STAT=stat)
      if(stat/=0) call cq_abort("Deallocation error in PulayMixSC: ",maxngrid, maxpulaySC)
      call reg_dealloc_mem(area_SC, 4*maxngrid*(1+maxpulaySC)+maxpulaySC*(maxpulaySC+1), type_dbl)
      return
    end subroutine dealloc_PulayMixSC

  end subroutine PulayMixSC

  subroutine PulayMixSCA(done,ndone,self_tol, reset_L, fixed_potential, vary_mu, n_L_iterations, &
       number_of_bands, L_tol, mu, total_energy,rho, size)
    use datatypes
    use numbers
    use PosTan
    use Pulay, ONLY: DoPulay2D
    use GenBlas
    use dimens, ONLY: n_my_grid_points, grid_point_volume
    use EarlySCMod, ONLY: get_new_rho, mixtwo
    use GenComms, ONLY: gsum, cq_abort, inode, ionode
    use io_module, ONLY: dump_charge, dump_charge2
    use hartree_module, ONLY: kerker
    use global_module, ONLY: ne_in_cell, iprint_SC, area_SC, flag_continue_on_SC_fail, flag_SCconverged
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use maxima_module, ONLY: maxngrid

    implicit none

    ! Passed variables
    logical :: done
    logical :: vary_mu, fixed_potential, reset_L

    integer :: size
    integer :: ndone
    integer :: n_L_iterations

    real(double) :: self_tol
    real(double) :: number_of_bands, L_tol, mixing, mu
    real(double) :: total_energy
    real(double), dimension(size) :: rho

    ! Local variables
    integer :: n_iters, n_pulay, npmod, pul_mx, i, j, ii, m, stat, thispulaySC
    real(double) :: R0, R1, E0, E1, dE, stLtol, tmp, max_neg, minA, maxA

    real(double), allocatable, dimension(:,:) :: rho_pul, R
    real(double), allocatable, dimension(:) :: rho1, resid, Kresid
    real(double), allocatable, dimension(:,:) :: Aij
    real(double), allocatable, dimension(:) :: alph

    character(len=11) :: digitstr = "01234567890"

!Reset Pulay Iterations  -- introduced by TM, Nov2007
    integer, parameter :: mx_fail = 3
    integer :: IterPulayReset=1, icounter_fail=0
    logical :: reset_Pulay=.false.
    real(double) :: R0_old

    allocate(rho_pul(maxngrid,maxpulaySC), R(maxngrid,maxpulaySC),rho1(maxngrid), &
         resid(maxngrid), Kresid(maxngrid), Aij(maxpulaySC, maxpulaySC), alph(maxpulaySC), STAT=stat)
    if(stat/=0) call cq_abort("Allocation error in PulayMixSC: ",maxngrid, maxpulaySC)
    call reg_alloc_mem(area_SC, maxngrid+2*maxngrid*(1+maxpulaySC)+maxpulaySC*(maxpulaySC+1), type_dbl)
    rho_pul = zero
    R = zero
    rho1 = zero
    resid = zero
    Kresid = zero
    Aij = zero
    alph = zero
    if(inode==ionode) write(io_lun,fmt='(8x,"Starting Pulay mixing, A = ",f6.3," q0= ",f7.4)') A, q0
    done = .false.
    n_iters = ndone
    m=1
    if(n_iters>=maxitersSC) then
       if(.NOT.flag_continue_on_SC_fail) call cq_abort('lateSC: too many iterations: ',n_iters, maxitersSC)
       flag_SCconverged = .false.
       done = .true.
       return
    endif
    ! m=1
    max_neg = 0.0_double
    do j=1,n_my_grid_points
       rho_pul(j,1) = rho(j)
       if(rho(j)<zero) max_neg = max(max_neg,abs(rho(j)))
    end do
    if(max_neg>0.0_double.AND.iprint_SC>3) &
         write(io_lun,fmt='(8x,"Init Max negative dens on node ",i5," is ",f8.3)') inode,max_neg
    call get_new_rho(.false., reset_L, fixed_potential, vary_mu, n_L_iterations, number_of_bands, L_tol, mu, &
         total_energy, rho, rho1, size)
    ! Evaluate residual
    do j=1,n_my_grid_points
       resid(j) = rho1(j) - rho(j)
    end do
    !tmp = grid_point_volume*asum(n_my_grid_points,resid,1)
    !write(io_lun,*) 'Sum of resid: ',tmp
    R0 = dot(n_my_grid_points,resid,1,resid,1)
    call gsum(R0)
    R0 = sqrt(grid_point_volume*R0)/ne_in_cell
    if(R0<self_tol) then
       if(inode==ionode) write(io_lun,'(8x,"Reached self-consistency tolerance")')
       done = .true.
       call dealloc_PulayMiXSCA
       return
    end if
    if(R0<EndLinearMixing) then
       if(inode==ionode) write(io_lun,'(8x,"Reached transition to LateSC")')
       call dealloc_PulayMiXSCA
       return
    end if
    !if(inode==ionode.AND.iprint_SC>=0) write(io_lun,*) 'Initial residual is ',R0
    if(inode==ionode) write(io_lun,fmt='(8x,"Pulay iteration ",i5," Residual is ",e12.5,/)') m,R0
    do j=1,n_my_grid_points
       R(j,1) = resid(j)
    end do
    ! Create new input charge
    Kresid = 0.0_double
    do j=1,n_my_grid_points
       Kresid(j) = resid(j)
    end do
    call kerker(Kresid,maxngrid,q0)
    !tmp = grid_point_volume*asum(n_my_grid_points,Kresid,1)
    !write(io_lun,*) 'Sum of Kresid: ',tmp
    max_neg = 0.0_double
    rho1 = zero
    rho1 = rho + Kresid
    !tmp = grid_point_volume*asum(n_my_grid_points,rho1,1)
    !call gsum(tmp)
    !write(io_lun,*) 'Sum of rho1: ',tmp
    !rho1 = ne_in_cell*rho1/tmp
    !write(io_lun,*) 'Calling mixtwo'
    call mixtwo(n_my_grid_points, .true., A, rho, rho1, resid)
    !tmp = grid_point_volume*asum(n_my_grid_points,resid,1)
    !write(io_lun,*) 'Sum of rho1: ',tmp
    rho = resid!*ne_in_cell/tmp
    !rho = rho1
    do j=1,n_my_grid_points
       if(rho(j)<0.0_double) then
          if(-rho(j)>max_neg) max_neg = -rho(j)
       !   rho(j) = 0.0_double
       end if
    end do
    if(max_neg>0.0_double.AND.iprint_SC>1) write(io_lun,*) 'First Max negative dens on node ',inode,max_neg
    ! Store deltarho
    n_iters = n_iters+1
 !Reset Pulay Iterations  -- introduced by TM, Nov2007
      R0_old = R0
      IterPulayReset=1
      icounter_fail = 0
    do m=2,maxitersSC
       ! calculate i (cyclical index for storing history) and pul_mx
       !ORI i = mod(m, maxpulaySC)
       !ORI if(i==0) i=maxpulaySC
       !ORI pul_mx = min(m, maxpulaySC)
      !Reset Version   TM Nov2007
       i = mod(m-IterPulayReset+1, maxpulaySC)
       if(i==0) i=maxpulaySC
       pul_mx = min(m-IterPulayReset+1, maxpulaySC)
       !if(inode==ionode) write(io_lun,fmt='(8x,"Pulay iteration ",i4)') m
       do j=1,n_my_grid_points
          rho_pul(j,i) = rho(j)
       end do
       ! Generate new charge and find residual
       call get_new_rho(.false., reset_L, fixed_potential, vary_mu, n_L_iterations, number_of_bands, L_tol, mu, &
            total_energy, rho, rho1, size)
       if(iprint_SC>1) call dump_charge(rho,size,inode)
       !call dump_charge2('ch'//digitstr(m:m),rho,size,inode)
       do j=1,n_my_grid_points
          resid(j) = rho1(j) - rho(j)
       end do
       !tmp = grid_point_volume*asum(n_my_grid_points,resid,1)
       !write(io_lun,*) 'Sum of resid: ',tmp
       R0 = dot(n_my_grid_points,resid,1,resid,1)
       call gsum(R0)
       R0 = sqrt(grid_point_volume*R0)/ne_in_cell
       if(inode==ionode) write(io_lun,fmt='(8x,"Pulay iteration ",i5," Residual is ",e12.5,/)') m,R0
       !call dump_charge2('re'//digitstr(m:m),resid,size,inode)
       !call dump_locps(resid,size,inode)
    !Reset Pulay Iterations
       if(R0 > R0_old ) icounter_fail = icounter_fail+1
       if(icounter_fail > mx_fail) then 
        if(inode == ionode) write(io_lun,*) ' Pulay iteration is reset !!  at ',m,'  th iteration'
        reset_Pulay = .true.
       endif 
        R0_old = R0
    !Reset Pulay Iterations  -- introduced by TM, Nov2007
       do j=1,n_my_grid_points
          R(j,i) = resid(j)
       end do
       if(R0<self_tol) then
          if(inode==ionode) write(io_lun,'(8x,"Reached self-consistent tolerance")')
          done = .true.
          call dealloc_PulayMiXSCA
          return
       end if
       if(R0<EndLinearMixing) then
          if(inode==ionode) write(io_lun,'(8x,"Reached transition to LateSC")')
          call dealloc_PulayMixSCA
          return
       end if
       ! now build new rho
       Aij = zero
       do ii=1,pul_mx
          R1 = dot(n_my_grid_points, R(:,ii), 1, R(:,ii),1)
          call gsum(R1)
          Aij(ii,ii) = R1
          if(ii>1) then
             do j=1,ii-1
                R1 = dot(n_my_grid_points, R(:,ii), 1, R(:,j),1)
                call gsum(R1)
                Aij(j,ii) = R1
                Aij(ii,j) = R1
             enddo
          end if
          !if(iprint_SC>2) write(io_lun,fmt='(5f16.12)') Aij(:,ii)
       enddo
       ! Solve to get alphas
       call DoPulay2D(Aij,alph,pul_mx,maxpulaySC,inode,ionode)
       if(inode==ionode.AND.iprint_SC>2) write(io_lun,*) 'alph: ',alph(1:pul_mx)
       ! Create delta rho - add on new rho
       rho = 0.0_double
       max_neg = 0.0_double
       do ii=1,pul_mx
          Kresid = zero
          Kresid(1:n_my_grid_points) = R(1:n_my_grid_points,ii)
          call kerker(Kresid,maxngrid,q0)
         !if(ii==pul_mx) then
          rho(1:n_my_grid_points) = rho(1:n_my_grid_points) + &
               alph(ii)*(rho_pul(1:n_my_grid_points,ii) + A*Kresid(1:n_my_grid_points))
         !else
         ! rho(1:n_my_grid_points) = rho(1:n_my_grid_points) + &
         !      alph(ii)*rho_pul(1:n_my_grid_points,ii) 
         !endif
       end do
       do j=1,n_my_grid_points
          if(rho(j)<0.0_double) then
             if(-rho(j)>max_neg) max_neg = -rho(j)
          end if
       end do
       if(max_neg>0.0_double.AND.iprint_SC>2) write(io_lun,*) i,' premix Max negative dens on node ',inode,max_neg
       n_iters = n_iters+1
     !Reset Pulay Iterations  -- introduced by TM, Nov2007
      if(reset_Pulay) then
        rho_pul(1:n_my_grid_points, 1) = rho_pul(1:n_my_grid_points, i)
        R(1:n_my_grid_points, 1) = R(1:n_my_grid_points, i)
        IterPulayReset = m
        icounter_fail=0
      endif
      reset_Pulay = .false.
     !Reset Pulay Iterations  -- introduced by TM, Nov2007
    end do
    ndone = n_iters
    call dealloc_PulayMixSCA
    return
    contains

      subroutine dealloc_PulayMixSCA

        implicit none

        integer :: stat

        deallocate(rho_pul, R,rho1, resid, Kresid, Aij, alph, STAT=stat)
        if(stat/=0) call cq_abort("Deallocation error in PulayMixSCA: ",maxngrid, maxpulaySC)
        call reg_dealloc_mem(area_SC, maxngrid+2*maxngrid*(1+maxpulaySC)+maxpulaySC*(maxpulaySC+1), type_dbl)
        return
      end subroutine dealloc_PulayMixSCA

  end subroutine PulayMixSCA

  subroutine PulayMixSCB(done,ndone,self_tol, reset_L, fixed_potential, vary_mu, n_L_iterations, &
       number_of_bands, L_tol, mu, total_energy,rho, size)
    use datatypes
    use numbers
    use PosTan
    use Pulay, ONLY: DoPulay2D
    use GenBlas
    use dimens, ONLY: n_my_grid_points, grid_point_volume
    use EarlySCMod, ONLY: get_new_rho, mixtwo
    use GenComms, ONLY: gsum, cq_abort, inode, ionode
    use io_module, ONLY: dump_charge, dump_charge2
    use hartree_module, ONLY: kerker
    use global_module, ONLY: ne_in_cell, flag_continue_on_SC_fail, flag_SCconverged
    use maxima_module, ONLY: maxngrid

    implicit none

    ! Passed variables
    logical :: done
    logical :: vary_mu, fixed_potential, reset_L

    integer :: size
    integer :: ndone
    integer :: n_L_iterations

    real(double) :: self_tol
    real(double) :: number_of_bands, L_tol, mixing, mu
    real(double) :: total_energy
    real(double), dimension(size) :: rho

    ! Local variables
    integer :: n_iters, n_pulay, npmod, pul_mx, i, j, ii, m, stat, thispulaySC
    real(double) :: R0, R1, E0, E1, dE, stLtol, tmp, max_neg

    real(double), allocatable, dimension(:,:) :: rho_pul, R
    real(double), allocatable, dimension(:) :: rho1, resid, Kresid
    real(double), allocatable, dimension(:,:) :: Aij
    real(double), allocatable, dimension(:) :: alph

    character(len=11) :: digitstr = "01234567890"

    allocate(rho_pul(maxngrid,maxpulaySC), R(maxngrid,maxpulaySC),rho1(maxngrid), &
         resid(maxngrid), Kresid(maxngrid), Aij(maxpulaySC, maxpulaySC), alph(maxpulaySC), STAT=stat)
    if(stat/=0) call cq_abort("Allocation error in PulayMixSC: ",maxngrid, maxpulaySC)
    rho_pul = zero
    R = zero
    rho1 = zero
    resid = zero
    Kresid = zero
    Aij = zero
    alph = zero
    if(inode==ionode) write(io_lun,fmt='(8x,"Starting Pulay mixing, A = ",f6.3," q0= ",f7.4)') A, q0
    done = .false.
    n_iters = ndone
    if(n_iters>=maxitersSC) then
       if(.NOT.flag_continue_on_SC_fail) call cq_abort('lateSC: too many iterations: ',n_iters, maxitersSC)
       flag_SCconverged = .false.
       done = .true.
       return
    endif
    ! m=1
    max_neg = 0.0_double
    do j=1,n_my_grid_points
       rho_pul(j,1) = rho(j)
    end do
    if(max_neg>0.0_double) write(io_lun,*) 'Init Max negative dens on node ',inode,max_neg
    call get_new_rho(.false., reset_L, fixed_potential, vary_mu, n_L_iterations, number_of_bands, L_tol, mu, &
         total_energy, rho, rho1, size)
    ! Evaluate residual
    do j=1,n_my_grid_points
       resid(j) = rho1(j) - rho(j)
    end do
    call kerker(resid,maxngrid,q0)
    !tmp = grid_point_volume*rsum(n_my_grid_points,resid,1)
    !write(io_lun,*) 'Sum of resid: ',tmp
    R0 = dot(n_my_grid_points,resid,1,resid,1)
    call gsum(R0)
    R0 = sqrt(grid_point_volume*R0)/ne_in_cell
    if(R0<self_tol) then
       if(inode==ionode) write(io_lun,'(8x,"Reached tolerance")')
       done = .true.
       call dealloc_PulayMixSCB
       return
    end if
    if(R0<EndLinearMixing) then
       if(inode==ionode) write(io_lun,'(8x,"Reached transition to LateSC")')
       call dealloc_PulayMixSCB
       return
    end if
    if(inode==ionode) write(io_lun,*) 'Residual is ',R0
    do j=1,n_my_grid_points
       R(j,1) = resid(j)
    end do
    rho1 = zero
    rho1 = rho + A*resid
    !tmp = grid_point_volume*rsum(n_my_grid_points,rho1,1)
    !call gsum(tmp)
    !write(io_lun,*) 'Sum of rho1: ',tmp
    !rho1 = ne_in_cell*rho1/tmp
    rho = rho1
    ! Store deltarho
    n_iters = n_iters+1
    m = 1
    i = 1
    do while (m<=maxitersSC)
       m = m+1
       i = i+1
       ! calculate i (cyclical index for storing history) and pul_mx
       i = mod(i, maxpulaySC)
       if(i==0) i=maxpulaySC
       ! Need something here to replace WORST residual
       pul_mx = min(m, maxpulaySC)
       if(inode==ionode) write(io_lun,fmt='(8x,"Pulay iteration ",i4," counters are ",2i4)') m, i, pul_mx
       do j=1,n_my_grid_points
          rho_pul(j,i) = rho(j)
       end do
       ! Generate new charge and find residual
       call get_new_rho(.false., reset_L, fixed_potential, vary_mu, n_L_iterations, number_of_bands, L_tol, mu, &
            total_energy, rho, rho1, size)
       !call dump_charge2('ch'//digitstr(m:m),rho,size,inode)
       do j=1,n_my_grid_points
          resid(j) = rho1(j) - rho(j)
       end do
       call kerker(resid,maxngrid,q0)
       R0 = dot(n_my_grid_points,resid,1,resid,1)
       call gsum(R0)
       R0 = sqrt(grid_point_volume*R0)/ne_in_cell
       if(inode==ionode) write(io_lun,*) 'Residual is ',R0
       !call dump_charge2('re'//digitstr(m:m),resid,size,inode)
       !call dump_locps(resid,size,inode)
       do j=1,n_my_grid_points
          R(j,i) = resid(j)
       end do
       if(R0<self_tol) then
          if(inode==ionode) write(io_lun,'(8x,"Reached tolerance")')
          done = .true.
          call dealloc_PulayMixSCB
          return
       end if
       if(R0<EndLinearMixing) then
          if(inode==ionode) write(io_lun,'(8x,"Reached transition to LateSC")')
          call dealloc_PulayMixSCB
          return
       end if
       ! now build new rho
       Aij = zero
       write(io_lun,*) 'A is :'
       do ii=1,pul_mx
          do j=1,pul_mx
             R1 = dot(n_my_grid_points, R(:,ii), 1, R(:,j),1)
             call gsum(R1)
             Aij(j,ii) = R1
          enddo
          write(io_lun,fmt='(5f16.12)') Aij(:,ii)
       enddo
       !if(minA/maxA<0.001_double) then
       !   write(io_lun,*) '
       !end if
       ! Solve to get alphas
       call DoPulay2D(Aij,alph,pul_mx,maxpulaySC,inode,ionode)
       if(inode==ionode)write(io_lun,*) 'alph: ',alph(1:pul_mx)
       ! Create delta rho - add on new rho
       rho = 0.0_double
       max_neg = 0.0_double
       do ii=1,pul_mx
          !Kresid = zero
          !Kresid(1:n_my_grid_points) = R(1:n_my_grid_points,ii)
          !call kerker(Kresid,maxngrid,q0)
          rho(1:n_my_grid_points) = rho(1:n_my_grid_points) + &
               alph(ii)*(rho_pul(1:n_my_grid_points,ii) + A*R(1:n_my_grid_points,ii))
               !alph(ii)*(rho_pul(1:n_my_grid_points,ii) + A*Kresid(1:n_my_grid_points))
       end do
       tmp = grid_point_volume*rsum(n_my_grid_points,rho,1)
       call gsum(tmp)
       write(io_lun,*) 'Sum of rho1: ',tmp
       rho = ne_in_cell*rho/tmp
       !rho = rho*ne_in_cell/tmp
       !do j=1,n_my_grid_points
          !if(rho(j)<0.0_double) then
          !   if(-rho(j)>max_neg) max_neg = -rho(j)
          !   rho(j) = 0.0_double
          !end if
       !end do
       if(max_neg>0.0_double) write(io_lun,*) i,' premix Max negative dens on node ',inode,max_neg
       n_iters = n_iters+1
    end do
    ndone = n_iters
    call dealloc_PulayMixSCB
    return
  contains

    subroutine dealloc_PulayMixSCB

      use GenComms, ONLY: cq_abort
      implicit none
      integer :: stat

      deallocate(rho_pul, R,rho1, resid, Kresid, Aij, alph, STAT=stat)
      if(stat/=0) call cq_abort("Deallocation error in PulayMixSCB: ",maxngrid, maxpulaySC)
      return
    end subroutine dealloc_PulayMixSCB

  end subroutine PulayMixSCB


end module SelfCon
