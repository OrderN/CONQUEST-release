! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module SelfCon_module
! ------------------------------------------------------------------------------
! Code area 5: self-consistency
! ------------------------------------------------------------------------------

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
!!   2015/06/08 lat
!!    - Added experimental backtrace
!!   2017/05/22 dave
!!    - Added minimum number of SCF iterations variable, minitersSC
!!   2019/08/16 14:01 dave
!!    Tidying: removing LinearMixSC and earlySC
!!   2019/08/13 16:46 jtlp with dave
!!    Implemented new residuals
!!   2019/10/24 11:52 dave
!!    Changed function calls to FindMinDM
!!   2019/12/30 tsuyoshi
!!    Introduced n_dumpSCF for dumping Kmatrix2.i99.p*****
!!  SOURCE
!!
module SelfCon

  use datatypes
  use global_module,          only: io_lun, area_SC
  use timer_module,           only: start_timer, stop_timer, cq_timer
  use timer_module,           only: start_backtrace, stop_backtrace
  use timer_stdclocks_module, only: tmr_std_chargescf, tmr_std_allocation

  implicit none

  ! Area identification
  integer, parameter, private :: area = 5

  ! These should all be read or dynamically allocated
  real(double), parameter :: thresh        = 2.0_double
  real(double), parameter :: InitialLambda = 1.0_double
  real(double), parameter :: ReduceLimit   = 0.5_double
  real(double), parameter :: crit_lin      = 0.1_double
  integer :: maxitersSC
  integer :: minitersSC
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
  ! Flags for new (2019/08) residual definitions
  logical, save :: flag_newresidual
  logical, save :: flag_newresid_abs


  real(double), dimension(2), save :: A ! A(spin) mixing factor for spin
  real(double), save :: q0
  real(double), save :: q1
  logical,      save :: flag_linear_mixing
  real(double), save :: EndLinearMixing

  ! Monitor change in energy from step to step of SCF cycle (for comparison with other minimisations)
  real(double), save :: dE_SCF

  integer,      save :: n_dumpSCF
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
  !!   2016/06/24 07:27 dave
  !!    Bug fix: removed comments hiding non-self-consistent exit
  !!  SOURCE
  !!
  subroutine new_SC_potl(record, self_tol, reset_L, fixed_potential, &
                         vary_mu, n_L_iterations, L_tol, total_energy, level)

    use datatypes
    use PosTan,          only: PulayC, PulayBeta, SCC, SCBeta, pos_tan, &
                               max_iters, SCE, SCR, fit_coeff
    use numbers
    use global_module,   only: iprint_SC, flag_self_consistent,         &
                               flag_SCconverged, IPRINT_TIME_THRES2,    &
                               nspin, spin_factor, min_layer, flag_mix_L_SC_min
    use H_matrix_module, only: get_H_matrix, get_output_energies
    use DMMin,           only: FindMinDM, dE_DMM
    use energy,          only: get_energy
    use GenComms,        only: inode, ionode, cq_abort
    use dimens,          only: n_my_grid_points
    use density_module,  only: density, get_electronic_density
    use maxima_module,   only: maxngrid
    use timer_module
    use io_module,       only: return_prefix
    use functions_on_grid, only: atomfns, H_on_atomfns

    implicit none

    ! Passed variables
    integer, optional :: level
    logical      :: record   ! Flags whether to record dE vs R
    logical      :: vary_mu, fixed_potential, reset_L
    integer      :: n_L_iterations
    real(double) :: self_tol
    real(double) :: L_tol
    real(double) :: total_energy


    ! Local variables
    real(double), dimension(nspin) :: electrons
    integer        :: ndone, i, nkeep, ndelta, stat, spin
    logical        :: done
    real(double)   :: SC_tol, DMM_tol, LastE
    type(cq_timer) :: tmr_l_tmp1
    type(cq_timer) :: backtrace_timer
    integer        :: backtrace_level
    character(len=12) :: subname = "new_SC_potl: "
    character(len=120) :: prefix

    prefix = return_prefix(subname, min_layer)
!****lat<$
    if (       present(level) ) backtrace_level = level+1
    if ( .not. present(level) ) backtrace_level = -10
    call start_backtrace(t=backtrace_timer,who='new_SC_potl',where=area,&
         level=backtrace_level,echo=.true.)
!****lat>$
    call start_timer(tmr_std_chargescf)

    ! Build H matrix *with NL and KE*
    call get_H_matrix(.true., fixed_potential, electrons, density, &
                      maxngrid, backtrace_level)

    ! The H matrix build is already timed on its own, so I leave out
    ! of the SC preliminaries (open to discussion)
    call start_timer(tmr_l_tmp1,WITH_LEVEL)
    if (.not. flag_self_consistent) then
       call stop_timer(tmr_std_chargescf)
       !min_layer = min_layer - 1
       call FindMinDM(n_L_iterations, vary_mu, L_tol, &
            reset_L, .false.)
       !if (.not.flag_mix_L_SC_min) then
       !call get_electronic_density(density, electrons, atomfns, &
       !     H_on_atomfns(1), inode, ionode, maxngrid, backtrace_level)
       !end if
       dE_SCF = dE_DMM
       call start_timer(tmr_std_chargescf)
       call get_output_energies(density, maxngrid)
       call get_energy(total_energy)
       call stop_print_timer(tmr_l_tmp1, "new_SC_potl (except H build)", &
                             IPRINT_TIME_THRES2)
       call stop_timer(tmr_std_chargescf)
       !min_layer = min_layer + 1
       return
    end if
    if (inode == ionode .and. iprint_SC + min_layer > 1) &
         write (io_lun, fmt='(4x,a,e12.5,/)') &
         trim(prefix)//" Starting.  Tolerance: ",self_tol
    if (record) then
       if (inode == ionode .and. iprint_SC + min_layer > 1) &
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
    flag_SCconverged = .true.
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
          call PulayMixSC_spin(done, ndone, SC_tol, reset_L,          &
                               fixed_potential, vary_mu,              &
                               n_L_iterations, DMM_tol, total_energy, &
                               density, maxngrid, backtrace_level)

       end if ! (flag_linear_mixing)

       call stop_print_timer(tmr_l_tmp1, "SCF iteration - pulay", &
                             IPRINT_TIME_THRES2)

       if (.not. done) then ! Late stage strategy
          call start_timer(tmr_l_tmp1, WITH_LEVEL)
          reset_L = .true. !.false.
          if (inode == ionode .and. iprint_SC + min_layer > 1) &
               write (io_lun, fmt='(10x,a)') '********** LateSC **********'

          call lateSC(record, done, ndone, SC_tol, reset_L,     &
                      fixed_potential, vary_mu, n_L_iterations, &
                      DMM_tol, total_energy, density, maxngrid, &
                      backtrace_level)

          call stop_print_timer(tmr_l_tmp1, "SCF iteration - Late", &
                                IPRINT_TIME_THRES2)
       end if
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
       if (inode == ionode .and. iprint_SC + min_layer > 1) then
          write (io_lun, *) '  List of residuals and energies'
          if (ndone > 0) then
             do i = 1, ndone
                write (io_lun, 7) i, SCR(i), SCE(i)
             end do
          end if
       end if
       call fit_coeff(SCC, SCBeta, SCE, SCR, ndone)
       if (inode == ionode .and. iprint_SC + min_layer > 1) &
            write (io_lun, 6) SCC, SCBeta
    end if

    call stop_print_timer(tmr_l_tmp1, "finishing SCF (fitting coefficients)", &
                          IPRINT_TIME_THRES2)
    call stop_timer(tmr_std_chargescf)

!****lat<$
    call stop_backtrace(t=backtrace_timer,who='new_SC_potl',echo=.true.)
!****lat>$

    return

6   format(8x,'dE to dR parameters - C: ',f15.8,' beta: ',f15.8)
7   format(8x,i4,2f15.8)

  end subroutine new_SC_potl
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
  !!   2015/06/08 lat
  !!   - Added experimental backtrace
  !!   2019/08/16 14:10 dave
  !!    get_new_rho is now in this module
  !!   2021/07/19 15:50 dave
  !!    Added dump charge at convergence (and at every iteration if iprint
  !!    is high enough)
  !!  SOURCE
  !!
  subroutine lateSC(record, done, ndone, self_tol, reset_L,          &
                    fixed_potential, vary_mu, n_L_iterations, L_tol, &
                    total_energy, rho, size, level)

    use datatypes
    use numbers
    use PosTan
    use Pulay,          only: DoPulay
    use GenBlas,        only: dot, rsum
    use dimens,         only: n_my_grid_points, grid_point_volume
    use GenComms,       only: gsum, cq_abort, inode, ionode
    use hartree_module, only: kerker
    use global_module,  only: ne_in_cell, area_SC,      &
                              flag_continue_on_SC_fail, &
                              flag_SCconverged,         &
                              flag_fix_spin_population, &
                              ne_spin_in_cell, nspin, spin_factor, iprint_SC
    use maxima_module,  only: maxngrid
    use memory_module,  only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use density_module, only: flag_DumpChargeDensity
    use io_module,      only: dump_charge

    implicit none

    ! Passed variables
    integer, optional :: level
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
    real(double), dimension(:),   allocatable :: rho_tot
    real(double), dimension(maxpulaySC,maxpulaySC,nspin) :: Aij
    real(double), dimension(maxpulaySC,nspin) :: alph
    real(double), dimension(nspin) :: R, tmp

    type(cq_timer)    :: backtrace_timer
    integer           :: backtrace_level

!****lat<$
    if (       present(level) ) backtrace_level = level+1
    if ( .not. present(level) ) backtrace_level = -10
    call start_backtrace(t=backtrace_timer,who='lateSC',&
         where=area,level=backtrace_level,echo=.true.)
!****lat>$

    allocate(rho_pul(maxngrid,maxpulaySC,nspin),   &
             resid_pul(maxngrid,maxpulaySC,nspin), &
             rho1(maxngrid,nspin), STAT=stat)
    if (stat /= 0) &
         call cq_abort("late_SC: Error alloc mem: ", maxngrid, maxpulaySC)
    call reg_alloc_mem(area_SC, (2*maxpulaySC+1)*nspin*maxngrid, type_dbl)
    if (flag_DumpChargeDensity .and. nspin==1) allocate(rho_tot(maxngrid))

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
                     maxngrid, backtrace_level)
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
    if (inode == ionode) write (io_lun, fmt='(10x,a,f16.6)') 'Residual is ', R0

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
            write (io_lun, fmt='(10x,a,i5)') '********** Late iter ', n_iters

       n_pulay = n_pulay + 1

       ! print out charge
       if (flag_DumpChargeDensity .and. iprint_SC > 2) then
          if (nspin == 1) then
             rho_tot(:) = spin_factor * rho(:,1)
             call dump_charge(rho_tot, n_my_grid_points, inode, spin=0)
          else
             call dump_charge(rho(:,1), n_my_grid_points, inode, spin=1)
             call dump_charge(rho(:,2), n_my_grid_points, inode, spin=2)
          end if
       end if

       ! Storage for pulay charges/residuals
       npmod = mod(n_pulay, maxpulaySC) + 1
       pul_mx = min(n_pulay + 1, maxpulaySC)
       if (inode == ionode) &
            write (io_lun, fmt='(10x,a,2i5)') 'npmod, pul_mx: ', npmod, pul_mx

       ! For the present output, find the residual (i.e. R_n^\prime)
       call get_new_rho(.false., reset_L, fixed_potential, vary_mu, &
                        n_L_iterations, L_tol, total_energy, rho,   &
                        rho1, maxngrid, backtrace_level)

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
                           resid_pul(:,i,spin), 1)
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

       call DoPulay(npmod, Aij, alph, pul_mx, maxpulaySC)

       !do spin = 1, nspin
       !   if (inode == ionode) &
       !        write (io_lun, *) 'alph (spin=', spin, '): ', &
       !                          alph(1:pul_mx,spin)
       !end do

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

          if (inode == ionode) write (io_lun, fmt='(10x,a)') 'Generating rho exactly'
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
                           rho1, maxngrid, backtrace_level)

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
               write (io_lun, fmt='(10x,a)') 'Generating rho by interpolation'
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
          if (inode == ionode) write (io_lun, fmt='(10x,a)') 'PANIC ! Residual increase !'
       end if
       !
       R0 = R1
       !
       if (R0 < self_tol) then
          done = .true.
          if (inode == ionode) write (io_lun, fmt='(10x,a)') 'Done ! Self-consistent'
       end if
       if (inode == ionode) write (io_lun, fmt='(10x,a,f16.6)') 'Residual is ', R0
       if (record) then
          SCE(n_iters) = total_energy
          SCR(n_iters) = R0
       end if
    end do

    ! print out charge
    if (flag_DumpChargeDensity) then
       if (nspin == 1) then
          rho_tot(:) = spin_factor * rho(:,1)
          call dump_charge(rho_tot, n_my_grid_points, inode, spin=0)
          deallocate(rho_tot)
       else
          call dump_charge(rho(:,1), n_my_grid_points, inode, spin=1)
          call dump_charge(rho(:,2), n_my_grid_points, inode, spin=2)
       end if
    end if

    ndone = n_iters

    if (inode == ionode) &
         write(io_lun, fmt='(10x,a,i5,a,f16.6)') 'Finishing lateSC after ', ndone, &
                          'iterations with residual of ', R0

    deallocate(rho_pul, resid_pul, rho1, STAT=stat)
    if (stat /= 0) &
         call cq_abort("late_SC: Error dealloc mem")
    call reg_dealloc_mem(area_SC, (2*maxpulaySC+1)*nspin*maxngrid, type_dbl)

!****lat<$
    call stop_backtrace(t=backtrace_timer,who='lateSC',echo=.true.)
!****lat>$

    return
  end subroutine lateSC
  !!***

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
  !!   2015/06/08 lat
  !!  - Added experimental backtrace
  !!   2019/08/13 16:53 jtlp with dave
  !!    Adding new residual definitions
  !!   2019/09/11 09:52 dave
  !!    Bug fix to changes in residual calculation
  !!   2019/12/30 tsuyoshi
  !!    flag_DumpChargeDensity is introduced to control dump_charge
  !!   2021/07/19 15:48 dave
  !!    Changed so that charge density is only output at every iteration
  !!    given iprint_SC>2.  Always output at SCF if flag is set.
  !! SOURCE
  !!
  subroutine PulayMixSC_spin(done, ndone, self_tol, reset_L, &
                             fixed_potential, vary_mu, n_L_iterations,&
                             L_tol, total_energy, rho, size, level)
    use datatypes
    use numbers
    use PosTan
    use GenBlas
    use dimens,         only: n_my_grid_points, grid_point_volume
    use GenComms,       only: gsum, cq_abort, inode, ionode, my_barrier
    use io_module,      only: dump_charge, return_prefix
    use hartree_module, only: kerker, kerker_and_wdmetric, wdmetric
    use global_module,  only: ne_in_cell, iprint_SC, area_SC,  &
                              flag_continue_on_SC_fail,        &
                              flag_SCconverged,                &
                              flag_fix_spin_population, nspin, &
                              spin_factor, min_layer
    use memory_module,  only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use maxima_module,  only: maxngrid
    use store_matrix,   only: dump_pos_and_matrices, unit_SCF_save
    use density_module, only: flag_DumpChargeDensity

    implicit none

    ! Passed variables
    integer,  optional :: level
    logical      :: done, vary_mu, fixed_potential, reset_L
    integer      :: size, ndone, n_L_iterations
    real(double) :: self_tol, L_tol, mixing, total_energy
    real(double), dimension(:,:) :: rho

    ! Local variables
    integer      :: n_iters, pul_mx, ii, iter, iPulay, stat
    real(double) :: R0, RA, RB, RC
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
    integer      :: spin, dotn
    real(double) :: R0_old, RA_old, RB_old, RC_old

    type(cq_timer)    :: backtrace_timer
    integer           :: backtrace_level
    character(len=12) :: subname = "PulayMixSC: "
    character(len=120) :: prefix

    prefix = return_prefix(subname, min_layer)
!****lat<$
    if (       present(level) ) backtrace_level = level+1
    if ( .not. present(level) ) backtrace_level = -10
    call start_backtrace(t=backtrace_timer,who='PulayMixSC_spin',&
         where=area,level=backtrace_level,echo=.true.)
!****lat>$

    ! allocate memories
    call allocate_PulayMixSC_spin

    ! initialise
    rho_pul = zero
    R_pul   = zero

    ! write out start information
    if (inode == ionode .and. iprint_SC + min_layer>0) then
       write (io_lun, '(4x,a,f6.3,a,f6.3)') &
             trim(prefix)//' Starting Pulay mixing, A_up = ', A(1), ' A_dn = ', A(2)
       if (nspin == 2) then
          if (flag_fix_spin_population) then
             write (io_lun, '(4x,a)') trim(prefix)//" Spin populations are fixed."
          else
             write (io_lun, '(4x,a)') trim(prefix)//" Spin populations are to be relaxed."
          end if
       else
          write (io_lun, '(4x,a)') trim(prefix)//" Spin non-polarised calculation."
       end if
       if (flag_Kerker) &
            write (io_lun, '(4x,a,f6.3)') &
                  trim(prefix)//' with Kerker preconditioning, q0 = ', q0
       if (flag_wdmetric) &
            write (io_lun, '(4x,a,f6.3)') &
                  trim(prefix)//' with wave dependent metric, q1 = ', q1
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

    dE_SCF = total_energy
    ! store rho and calculate residuals and store in pulay history slot 1
    call update_pulay_history(1, rho, reset_L, fixed_potential,     &
                              vary_mu, n_L_iterations, L_tol,       &
                              total_energy, rho_pul, R_pul, KR_pul, &
                              Rcov_pul, backtrace_level)
    dE_SCF = total_energy - dE_SCF

    ! Evaluate magnitute of residual, note do not include cross terms
    RA = zero
    do spin = 1, nspin
       RA = RA + spin_factor * &
            dot(n_my_grid_points, R_pul(:,1,spin), 1, R_pul(:,1,spin), 1)
    end do
    call gsum(RA)
    RA = sqrt(grid_point_volume * RA) / ne_in_cell
    ! New, absolute
    RB = zero
    do spin = 1, nspin
       RB = RB + spin_factor * sum(abs(R_pul(:,1,spin)))
    end do
    call gsum(RB)
    RB = grid_point_volume * RB
    ! New, relative
    RC = RB / ne_in_cell

    ! Set residual
    if ( .not. flag_newresidual) then
       R0 = RA
    else
       if(flag_newresid_abs) then
          R0 = RB
       else
          R0 = RC
       end if
    end if
    ! print residual information
    if (inode==ionode) then
       if (iprint_SC + min_layer>2) then
          write (io_lun, '(4x,a,i5,a,e12.5)') &
               trim(prefix)//' Pulay iteration ', iter, ' RMS residual:             ', RA
          write (io_lun, '(4x,a,i5,a,e12.5)') &
               trim(prefix)//' Pulay iteration ', iter, ' absolute residual (tot) : ', RB
          write (io_lun, '(4x,a,i5,a,e12.5)') &
               trim(prefix)//' Pulay iteration ', iter, ' absolute residual (frac): ', RC
      else if(iprint_SC + min_layer>0) then
          write (io_lun, '(4x,a,i5,a,e12.5)') &
               trim(prefix)//' Pulay iteration ', iter, ' residual:             ', R0
      end if
    end if
    ! check if they have reached tolerance
    if (R0 < self_tol .AND. iter >= minitersSC) then ! If we've done minimum number
       if (inode == ionode .and. iprint_SC + min_layer>=0) &
            write (io_lun,fmt='(4x,a,e12.5, &
            &" after ",i6," iterations")') trim(prefix)//" reached SCF residual of ", &
            R0, iter
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

       ! print out charge
       !  2019Dec30 tsuyoshi: flag_DumpChargeDensity is introduced, but ...
       !                      is it okay with iprint_SC + min_layer > 1 ?
       if (flag_DumpChargeDensity .and. iprint_SC + min_layer > 1) then
          if (nspin == 1) then
             rho_tot(:) = spin_factor * rho(:,1)
             call dump_charge(rho_tot, n_my_grid_points, inode, spin=0)
          else
             call dump_charge(rho(:,1), n_my_grid_points, inode, spin=1)
             call dump_charge(rho(:,2), n_my_grid_points, inode, spin=2)
          end if
       end if

       dE_SCF = total_energy
       ! calculate cyclic index for storing pulay history
       iPulay = mod(iter - IterPulayReset + 1, maxpulaySC)
       if (iPulay == 0) iPulay = maxpulaySC
       ! calculated the number of pulay histories stored
       pul_mx = min(iter - IterPulayReset + 1, maxpulaySC)

       ! Calculate residuals and update pulay history (store in iPulay-th slot)
       call update_pulay_history(iPulay, rho, reset_L, fixed_potential, &
                                 vary_mu, n_L_iterations, L_tol,        &
                                 total_energy, rho_pul, R_pul, KR_pul,  &
                                 Rcov_pul, backtrace_level)
       dE_SCF = total_energy - dE_SCF

       ! Evaluate magnitute of residual, note no cross terms
       RA = zero
       RB = zero
       RC = zero
       ! Original
       do spin = 1, nspin
          RA = RA + spin_factor * &
               dot(n_my_grid_points, R_pul(:,iPulay,spin), 1, &
               R_pul(:,iPulay,spin), 1)
       end do
       call gsum(RA)
       RA = sqrt(grid_point_volume * RA) / ne_in_cell
       ! New, absolute
       do spin = 1, nspin
          RB = RB + spin_factor * sum(abs(R_pul(:,iPulay,spin)))
       end do
       call gsum(RB)
       RB = grid_point_volume * RB
       ! New, relative
       RC = RB / ne_in_cell
       ! Set residual
       if ( .not. flag_newresidual) then
          R0 = RA
       else
          if(flag_newresid_abs) then
             R0 = RC
          else
             R0 = RB
          end if
       end if
       ! print residual information
       if (inode==ionode) then
          if (iprint_SC + min_layer>2) then
             write (io_lun, '(4x,a,i5,a,e12.5)') &
                  trim(prefix)//' Pulay iteration ', iter, ' RMS residual:             ', RA
             write (io_lun, '(4x,a,i5,a,e12.5)') &
                  trim(prefix)//' Pulay iteration ', iter, ' absolute residual (tot) : ', RB
             write (io_lun, '(4x,a,i5,a,e12.5)') &
                  trim(prefix)//' Pulay iteration ', iter, ' absolute residual (frac): ', RC
          else if(iprint_SC + min_layer>0) then
             write (io_lun, '(4x,a,i5,a,e12.5)') &
                  trim(prefix)//' Pulay iteration ', iter, ' residual:             ', R0
          end if
       end if
       ! check if they have reached tolerance
       if (R0 < self_tol .AND. iter >= minitersSC) then ! Passed minimum number of iterations
          if (inode == ionode .and. iprint_SC + min_layer>=0) &
               write (io_lun,fmt='(4x,a,e12.5, &
               &" after ",i6," iterations")') trim(prefix)//" Reached SCF tolerance of ", &
               R0, iter
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
          if (inode == ionode .and. iprint_SC + min_layer>0) &
               write (io_lun, fmt='(4x,a,i3,a)') &
               trim(prefix)//' Pulay iteration is reset !!  at ', iter, &
               ' th iteration'
          reset_Pulay = .true.
       end if
       R0_old = R0
       ! get optimum rho mixed from rho_pul, R_pul and Rcov_pul
       call get_pulay_optimal_rho(iPulay, rho, pul_mx, A, rho_pul, &
                                  R_pul, KR_pul, Rcov_pul, backtrace_level)

       ! 2019Dec30 tsuyoshi
       ! Dump Kmatrix every n_dumpSCF, if n_dumpSCF > 0
       if(n_dumpSCF > 0 .and. mod(n_iters,n_dumpSCF)==1) then
        call dump_pos_and_matrices(index=unit_SCF_save)
       endif

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
    call stop_backtrace(t=backtrace_timer,who='PulayMixSC_spin',echo=.true.)
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
  !!   2015/06/08 lat
  !!   - Added experimental backtrace
  !! SOURCE
  !!
  subroutine get_pulay_optimal_rho(iPulay, rho_opt, pul_max,       &
                                   mix_factor, rho_pul, resid_pul, &
                                   k_resid_pul, cov_resid_pul, level)
    use datatypes
    use GenBlas
    use global_module,  only: nspin, flag_fix_spin_population
    use GenComms,       only: gsum, inode, ionode
    use Pulay,          only: DoPulay
    use dimens,         only: n_my_grid_points
    use maxima_module,  only: maxngrid

    implicit none
 
    ! Passed parameters
    integer, optional   :: level    
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
    type(cq_timer) :: backtrace_timer
    integer        :: backtrace_level

!****lat<$
    if (       present(level) ) backtrace_level = level+1
    if ( .not. present(level) ) backtrace_level = -10
    call start_backtrace(t=backtrace_timer,who='get_pulay_optimal_rho',&
         where=area,level=backtrace_level,echo=.true.)
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
    call DoPulay(iPulay, Aij, alpha, pul_max, maxpulaySC)

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
    call stop_backtrace(t=backtrace_timer,who='get_pulay_optimal_rho',echo=.true.)
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
  !!   2015/06/08 lat
  !!    - Added experimental backtrace
  !!   2019/08/16 14:10 dave
  !!    get_new_rho is now in this module
  !! SOURCE
  !!
  subroutine update_pulay_history(iPulay, rho_in, reset_L,             &
                                  fixed_potential, vary_mu,            &
                                  n_L_iterations, L_tol, total_energy, &
                                  rho_pul, resid_pul, k_resid_pul,     &
                                  cov_resid_pul, level)
    use datatypes
    use GenBlas
    use dimens,         only: n_my_grid_points, grid_point_volume
    use GenComms,       only: gsum, cq_abort
    use hartree_module, only: kerker, kerker_and_wdmetric, wdmetric
    use global_module,  only: nspin, flag_fix_spin_population
    use maxima_module,  only: maxngrid
    use memory_module,  only: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    ! Passed parameters
    integer,      optional    :: level 
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
    type(cq_timer) :: backtrace_timer
    integer        :: backtrace_level

!****lat<$
    if (       present(level) ) backtrace_level = level+1
    if ( .not. present(level) ) backtrace_level = -10
    call start_backtrace(t=backtrace_timer,who='update_pulay_history',&
         where=area,level=backtrace_level,echo=.true.)
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
                     rho_out, maxngrid, backtrace_level)
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
    call stop_backtrace(t=backtrace_timer,who='update_pulay_history',echo=.true.)
!****lat>$

    return
  end subroutine update_pulay_history
  !*****

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
  !!   2020/08/24 11:21 dave
  !!    Add option to allow LFD at each SCF step
  !!   2021/07/26 11:49 dave
  !!    Fix module use for get_H_matrix
  !!   2021/07/28 10:13 dave
  !!    Correctly calculate DFT energy (Ha, XC and local contributions use output density)
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
                                 nspin, spin_factor, flag_diagonalisation, flag_LFD, flag_Multisite
    use H_matrix_module,   only: get_H_matrix, get_output_energies
    use S_matrix_module,   only: get_S_matrix
    !use DiagModule,        only: diagon
    use energy,            only: get_energy
    use functions_on_grid, only: atomfns, allocate_temp_fn_on_grid, &
                                 free_temp_fn_on_grid
    use density_module,    only: get_electronic_density
    use GenComms,          only: inode, ionode
    use cdft_module,       only: cdft_min
    use multisiteSF_module, only: initial_SFcoeff, flag_mix_LFD_SCF

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

    if(flag_Multisite .and. flag_LFD .and. flag_mix_LFD_SCF) then
       call initial_SFcoeff(.false.,.false.,fixed_potential,.false.)
       call get_S_matrix(inode, ionode, build_AtomF_matrix=.false.)
       call get_H_matrix(.false., fixed_potential, electrons, rhoin, &
            size, level=backtrace_level,build_AtomF_matrix=.false.)
    end if
    if (flag_perform_cDFT) then
       call cdft_min(reset_L, fixed_potential, vary_mu, &
                     n_CG_L_iterations, tolerance, total_energy)
    else
       Ltol = tolerance
       ! Find minimum density matrix
       call stop_timer(tmr_std_chargescf)
       call FindMinDM(n_CG_L_iterations, vary_mu, Ltol, &
                      reset_L, record, backtrace_level)
       call start_timer(tmr_std_chargescf)

       ! If we're using O(N), we only have L, and we need K
       ! if diagonalisation, we have K
       if (.not. flag_diagonalisation) then
          call LNV_matrix_multiply(electrons, energy_spin, doK,      &
                                   dontM1,  dontM2, dontM3, dontM4,  &
                                   dontphi, dontE, level=backtrace_level)
       end if
    end if ! if (flag_perform_cDFT) then

    ! And get the output density
    temp_supp_fn = allocate_temp_fn_on_grid(atomf)
    call stop_timer(tmr_std_chargescf) ! This routine is always call within area 5
    call get_electronic_density(rhoout, electrons, atomfns, &
                                temp_supp_fn, inode, ionode, size, backtrace_level)
    call start_timer(tmr_std_chargescf)
    call free_temp_fn_on_grid(temp_supp_fn)

    ! Now build DFT energy from output charge
    ! Find Hartree, XC and local PS (i.e. NA) energies with output density
    call get_output_energies(rhoout, size)
    call get_energy(total_energy, .true., backtrace_level) ! Output DFT energy

!****lat<$
    call stop_backtrace(t=backtrace_timer,who='get_new_rho',echo=.true.)
!****lat>$

    return

  end subroutine get_new_rho
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
  !!   2018/11/13 17:30 nakata
  !!   - Changed matS to be spin_SF dependent
  !!   2019/11/21 nakata 
  !!   - Introduced spin-up/down atomic charge
  !!   2019/11/21 10:51 dave
  !!   - Tweak to use gsum for 2D array for charge
  subroutine get_atomic_charge()

    use datatypes
    use numbers
    use global_module,  only: ni_in_cell, area_SC, nspin, spin_factor, &
                              flag_SpinDependentSF
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
    integer :: chun, stat, n,l, glob_ind, spin, spin_SF
    integer, dimension(nspin) :: temp_mat
    real(double), dimension(:,:), allocatable :: charge, node_charge

    call start_timer(tmr_std_chargescf)
    ! prepare arrays and suitable matrices
    l = bundle%n_prim
    call start_timer(tmr_std_allocation)
    allocate(charge(ni_in_cell,nspin), STAT=stat)
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
    spin_SF = 1
    do spin = 1, nspin
       if (flag_SpinDependentSF) spin_SF = spin
       temp_mat(spin) = allocate_temp_matrix(Srange, 0)
       call matrix_sum(zero, temp_mat(spin), one, matK(spin))
       call atom_trace(temp_mat(spin), matS(spin_SF), l, node_charge(:,spin))
       ! sum from the node_charge into the total charge array
       do n = 1, l
          glob_ind = atoms_on_node(n, inode)
          charge(glob_ind,spin) = charge(glob_ind,spin) + &
                                  spin_factor * node_charge(n,spin)
       end do
    end do
    call gsum(charge, ni_in_cell, nspin)

    ! output
    if (inode == ionode) then
       ! write (*, *) 'Writing charge on individual atoms...'
       call io_assign(chun)
       open(unit = chun, file='AtomCharge.dat')
       if (nspin.eq.1) then
          do n = 1, ni_in_cell
             write (chun, fmt='(f15.10)') charge(n,1)
          end do
       else
          do n = 1, ni_in_cell
             write (chun, fmt='(3f15.10)') charge(n,1)+charge(n,2), charge(n,1), charge(n,2)
          end do
       endif
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

