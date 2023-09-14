! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module DMMin
! ------------------------------------------------------------------------------
! Code area 4: Density matrix
! ------------------------------------------------------------------------------

!!****h* Conquest/DMMin
!!  NAME
!!   DMMin
!!  PURPOSE
!!   Contains all routines related to minimisation of energy with
!!   respect to elements of density matrix
!!
!!   Contains initialisation, early stage (steepest descents) and late
!!   stage (GR-Pulay - D.R.Bowler and M.J.Gillan, Chem. Phys.
!!   Lett. 325, 473 (2000)) techniques as well as a call to a
!!   diagonalisation routine for exact solution
!!  USES
!!   datatypes, DiagModule, GenBlas, GenComms, global_module,
!!   logicals, matrix_data, maxima_module, McWeeny, multiply_module,
!!   mult_module, numbers, PosTan, primary_module, Pulay
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   11/12/00
!!  MODIFICATION HISTORY
!!   05/06/2001 dave
!!    ROBODoc header, RCS Id and Log tags and changed all stops
!!    to cq_abort calls
!!   06/06/2001 dave
!!    Bug fixes - added temp_M to use matrix_data and closed
!!    brackets on cq_abort call
!!   08/06/2001 dave
!!    Changed all occurences of dgsum to gsum
!!   29/05/2002 dave
!!    Minor changes to shift H back after DM minimisation, and to
!!    switch between O(N) and diagonalisation
!!   14:23, 26/02/2003 drb
!!    Added electron number check to linear part of
!!    correct_electron_number
!!   16:43, 10/05/2005 dave
!!    Various small changes throughout code: mainly bug fixes
!!   2008/02/01 17:46 dave
!!    Changes to write output to file not stdout
!!   2008/05/22 ast
!!    Added timers
!!   2012/03/01 L.Tong
!!    Added interface for lineMinL
!!   2013/03/06 L.Tong
!!    Added module globals maxpulaystepDMM and minpulaystepDMM
!!    to allow user to control the pulay mixing stepsize bracket
!!    in lateDM
!!   2014/09/15 18:30 lat
!!    fixed call start/stop_timer to timer_module (not timer_stdlocks_module !)
!!   2019/10/24 08:30 dave
!!    Added inode and ionode as use-associated from GenComms in module header
!!    for efficiency and best practice; also iprint_DM
!!   2021/08/02 14:46 dave
!!    Added dE_DMM for comparison to structural optimisation
!!  SOURCE
!!
module DMMin

  use datatypes
  use global_module,          only: io_lun, area_DM, iprint_DM
  use timer_module,           only: start_timer, stop_timer, cq_timer
  use timer_module,           only: start_backtrace, stop_backtrace
  use timer_stdclocks_module, only: tmr_std_densitymat
  use GenComms,               ONLY: inode, ionode
  
  ! Area identification
  integer, parameter, private :: area = 4
  integer           :: maxpulayDMM
  integer,     save :: n_dumpL = 0     ! Redefined: 2019Dec30 tsuyoshi
  real(double)      :: LinTol_DMM
  real(double)      :: minpulaystepDMM ! min step size for pulay minimisation
  real(double)      :: maxpulaystepDMM ! max step size for pulay minimisation
  real(double), save:: dE_DMM
  !integer, parameter :: mx_pulay = 5
  !real(double), parameter :: LinTol = 0.1_double
  
!!***

contains

  ! -----------------------------------------------------------
  ! Subroutine FindMinDM
  ! -----------------------------------------------------------

  !!****f* DMMin/FindMinDM *
  !!
  !!  NAME
  !!   FindMinDM - control routine for energy minimisation wrt DM
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Controls the minimisation of energy with respect to the
  !!   density matrix.  Calls McWeeny for initialisation, and
  !!   early and late stage techniques for minimisation.
  !!
  !!   If the early stage finds an inflexion during exact line
  !!   minimisation, then the code immediately bounces back to
  !!   McWeeny initialisation.  If the residual increases during
  !!   GRP minimisation (which it shouldn't !) early stage is then
  !!   used.
  !!  INPUTS
  !!   logical :: resetL - determines whether it's necessary to use McWeeny
  !!   logical :: record - flag for recording E vs residual
  !!  USES
  !!   datatypes, numbers, global_module, maxima_module, matrix_data,
  !!   mult_module, multiply_module, McWeeny, PosTan
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   11/12/00
  !!  MODIFICATION HISTORY
  !!   05/06/2001 dave
  !!    Added ROBODoc header
  !!   20/06/2001 dave
  !!    Moved fitting of C and beta into PosTan and called fit_coeff
  !!   29/05/2002 dave
  !!    Added option to use diagonalisation to find DM
  !!   2004/10/28 drb
  !!    Added dump_matrix call at the end of minimisation
  !!   16:43, 10/05/2005 dave
  !!    Added local dataM3
  !!   2008/05/22 ast
  !!    Added timer
  !!   2011/07/22 L.Tong
  !!    Added spin polarisation
  !!   2011/09/13 L.Tong
  !!    Removed obsolete dependence on number_of_bands, and instead use
  !!    the global module variable ne_in_cell for the total number of
  !!    electrons
  !!   2012/03/08 L.Tong
  !!    - Major rewrite of spin implementation
  !!    - Removed redundant input parameter real(double) mu
  !!   2012/06/24 L.Tong
  !!    - Added control of whether to dump L using flag_dump_L
  !!   2013/08/20 M.Arita
  !!    - Changed calls for dumping L-matrix
  !!    - Added flag to start from EarlyDM
  !!   2013/12/03 M.Arita
  !!    - Added calls for writing out matrices for XL-BOMD
  !!    - Added call for writing out InfoGlobal.dat
  !!      (These used to be called at get_E_and_F)
  !!   2014/01/21 lat
  !!    - Added call to exx
  !!   2014/09/15 18:30 lat
  !!    - Fixed call start/stop_timer to timer_module (not timer_stdlocks_module !)
  !!   2015/06/08 lat
  !!    - Added experimental backtrace
  !!   2017/02/23 dave
  !!    - Changing location of diagon flag from DiagModule to global and name to flag_diagonalisation
  !!   2017/05/09 dave
  !!    Added code to dump L matrix for both spin channels
  !!   2017/05/11 dave
  !!    Added code to dump X matrix (and associated) for both spin channels
  !!   2017/11/06 dave
  !!    Added dump of K matrix
  !!   2018/01/22 tsuyoshi (with dave)
  !!    Initial changes for atom updates
  !!   2018/11/13 17:30 nakata
  !!    Changed matS, matT and matTtran to be spin_SF dependent
  !!   2019/05/27 tsuyoshi 
  !!    Debug for spin=2, propagateL (not propagateX)
  !!  SOURCE
  !!
  subroutine FindMinDM(n_L_iterations, vary_mu, tolerance, &
                       resetL, record, level)

    use datatypes
    use numbers
    use global_module, only: IPRINT_TIME_THRES1,             &
                             nspin, nspin_SF, flag_fix_spin_population, &
                             ne_in_cell, ne_spin_in_cell, flag_dump_L,  &
                             flag_SkipEarlyDM, flag_XLBOMD,             &
                             flag_propagateX, flag_dissipation,         &
                             integratorXL, runtype, flag_exx,           &
                             flag_diagonalisation, min_layer, flag_DM_converged
    use mult_module,   only: matrix_transpose, matT, matTtran, matL,    &
                             matS, matrix_sum, matK
    use McWeeny,       only: InitMcW, McWMin
    use PosTan,        only: max_iters, cscale, PulayE, PulayR, PulayC, &
                             PulayBeta, pos_tan, fit_coeff
    use DiagModule,    only: FindEvals
    use io_module,     only: dump_matrix, return_prefix
    use energy,        only: entropy
    use timer_module,  only: cq_timer, start_timer, stop_print_timer,   &
                             WITH_LEVEL
    use matrix_data,   only: Lrange, Srange, LSrange, Hrange
    !use XLBOMD_module, only: matX, matXvel, dump_XL
    use store_matrix,  only: dump_XL
    use mult_module,   only: matXL, matXLvel

    use exx_kernel_default, only: get_X_matrix

    implicit none

    ! Passed variables
    integer, optional :: level
    integer      :: n_L_iterations
    real(double) :: tolerance
    logical      :: resetL, vary_mu, record

    ! Local variables
    logical                 :: done, inflex, problem, early
    integer                 :: ndone, nkeep, ndelta, i
    real(double)            :: LastE
    real(double), parameter :: eprec = 1.0e-12_double
    real(double), save      :: delta_e
    type(cq_timer)          :: tmr_l_iter
    type(cq_timer)          :: backtrace_timer
    integer                 :: backtrace_level
    !TM 2010.Nov.06
    integer                 :: niter = 0
    integer, parameter      :: niter_max = 10
    character(len=12) :: subname = "FindDM: "
    character(len=120) :: prefix

!****lat<$
    if (       present(level) ) backtrace_level = level+1
    if ( .not. present(level) ) backtrace_level = -10
    call start_backtrace(t=backtrace_timer,who='FindMinDM',&
         where=area,level=backtrace_level,echo=.true.)
!****lat>$

    call start_timer(tmr_std_densitymat)
    prefix = return_prefix(subname, min_layer)

    entropy = zero

    if (flag_diagonalisation) then ! Use exact diagonalisation to get K
       call FindEvals(ne_spin_in_cell)
    else
       if (inode == ionode .and. iprint_DM + min_layer >0) then
          write(io_lun,'(/4x,a,f15.10,2x,a,i5)') &
               trim(prefix)//' tol: ', tolerance, &
               'number of L iterations: ', n_L_iterations
       end if

       problem = .false.
       inflex  = .false.
       done    = .false.
       early   = .true.  ! Do early DM iterations by default
       if (flag_SkipEarlyDM) then
          early = .false.
       end if

       ! Start with the transpose of S^-1
       call matrix_transpose(matT(1), matTtran(1))
       if (nspin_SF==2) call matrix_transpose(matT(2), matTtran(2))

       ! Now minimise the energy
       ndone = 0
       niter = 0
       do while (.not. done)
          call start_timer(tmr_l_iter, WITH_LEVEL)
          if (resetL .or. inflex) then ! Reset L to McW
             call InitMcW
             call McWMin(n_L_iterations, delta_e)
             early = .true.
             resetL = .false.    !2010.Nov.06 TM
          end if
          if (early .or. problem) then
             inflex = .false.
             call earlyDM(ndone, n_L_iterations, delta_e, done, vary_mu, &
                  tolerance, inflex, record)
          end if
          if ((.not. done) .and. (.not. inflex)) then ! Continue on to late stage
             call lateDM(ndone, n_L_iterations, done, delta_e, vary_mu, &
                  tolerance, record)
          end if
          niter = niter + 1
          if (problem) resetL = .true.
          !ORI problem = .true. ! If we're not done, then this kicks us back to early
          if (niter > niter_max) then
             problem = .true. ! If we're not done, then this kicks us back to early
             niter = 0
          end if
          call stop_print_timer(tmr_l_iter, "FindMinDM iteration", &
               IPRINT_TIME_THRES1)
       end do ! end of do while (.not. done)
    end if
    flag_DM_converged = .true.
    if (record) then
       if (inode == ionode .and. iprint_DM + min_layer > 2) then
          write(io_lun,*) '  List of residuals and energies'
          do i = 1, ndone
             write(io_lun, 7) i, PulayR(i), PulayE(i)
          end do
       end if
       call fit_coeff(PulayC, PulayBeta, PulayE, PulayR, ndone)
       if (inode == ionode) write (io_lun, 6) PulayC, PulayBeta
    end if

    call stop_timer(tmr_std_densitymat)

!****lat<$
    call stop_backtrace(t=backtrace_timer,who='FindMinDM',echo=.true.)
!****lat>$

    return

6   format(2x,'dE to dR parameters - C: ',f15.8,' beta: ',f15.8)
7   format(2x,i4,2f15.8)
8   format(2x,'Welcome to minimise_energy. Tolerance is ',f15.8)

  end subroutine FindMinDM
  !!***

  ! -----------------------------------------------------------
  ! Subroutine earlyDM
  ! -----------------------------------------------------------

  !!****f* DMMin/earlyDM *
  !!
  !!  NAME
  !!   earlyDM - early stage DM min
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Performs steepest descent minimisation of energy wrt density
  !!   matrix elements, until the change in gradient is linear in
  !!   density matrix (at which point CG or GR-Pulay - late stage in
  !!   other words) can be used
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   11/12/00
  !!  MODIFICATION HISTORY
  !!   05/06/2001 dave
  !!    Added ROBODoc header, removed calls to stop
  !!   08/06/2001 dave
  !!    Changed dgsum to gsum from GenComms
  !!   29/05/2002 dave
  !!    Added call to move H back (i.e. undo the H-mu.S shift) after
  !!    minimisation
  !!   2004/10/28 drb
  !!    Removed shift of H
  !!   10:09, 13/02/2006 drb
  !!    Removed all explicit references to data_ variables and rewrote
  !!    in terms of new matrix routines
  !!   2009/07/08 16:41 dave
  !!    Introduced atom-based tolerance (just divide by ni_in_cell)
  !!   2011/08/14 L.Tong
  !!    Added spin polarisation
  !!    Added matphi_dn
  !!   2011/08/25 L.Tong
  !!    Removed the redundant parameter number_of_bands
  !!   2012/03/12 L.Tong
  !!    - Major rewrite of spin implementation
  !!    - Removed redundant input parameter real(double) mu
  !!   2018/11/13 17:30 nakata
  !!    Changed matT to be spin_SF dependent
  !!  SOURCE
  !!
  subroutine earlyDM(ndone, n_L_iterations, delta_e, done, vary_mu, &
                     tolerance, inflex, record)

    use datatypes
    use logicals
    use numbers
    use matrix_data,    only: Lrange, TLrange
    use mult_module,    only: matT, matphi, mult, T_L_TL, TL_T_L,      &
                              matrix_product_trace,                    &
                              LNV_matrix_multiply,                     &
                              allocate_temp_matrix, free_temp_matrix,  &
                              matrix_sum, matrix_product
    use primary_module, only: bundle
    use PosTan,         only: PulayR, PulayE
    use GenComms,       only: cq_abort, gsum, cq_warn
    use global_module,  only: IPRINT_TIME_THRES1, min_layer,           &
                              ni_in_cell, flag_global_tolerance,       &
                              nspin, spin_factor,                      &
                              flag_fix_spin_population, flag_SpinDependentSF
    use timer_module,   only: cq_timer, start_timer, stop_print_timer, &
                              WITH_LEVEL
    use io_module,      only: return_prefix

    implicit none

    ! Passed variables
    real(double) :: tolerance, delta_e
    integer      :: n_L_iterations, ndone
    logical      :: inflex, vary_mu, done, record

    ! Local variables
    integer                        :: n_iter, length, spin, spin_SF
    integer,      dimension(nspin) :: matM3, matSM3, matSphi, mat_temp, mat_search
    real(double), dimension(nspin) :: energy0, energy1, electrons
    real(double)                   :: energy0_tot, energy1_tot, electrons_tot
    real(double), dimension(nspin) :: e_dot_n, n_dot_n
    real(double)                   :: e_dot_n_tot, n_dot_n_tot
    real(double)                   :: g0, g1
    real(double)                   :: interpG, zeta, direct_sum_factor
    type(cq_timer)                 :: tmr_l_tmp1, tmr_l_iter
    character(len=12) :: subname = "earlyDM: "
    character(len=120) :: prefix

    ! integer :: matM3, matSM3, matSphi, mat_temp, mat_search
    ! integer :: matM3_dn, matSM3_dn, mat_search_dn, matSphi_dn
    ! real(double) :: energy0, energy0_up, energy0_dn
    ! real(double) :: energy1, energy1_up, energy1_dn
    ! real(double) :: electrons, electrons_up, electrons_dn
    ! real(double) :: e_dot_n, e_dot_n_up, e_dot_n_dn
    ! real(double) :: n_dot_n, n_dot_n_up, n_dot_n_dn
    ! real(double) :: g0, g1
    ! real(double) ::  interpG,  zeta, direct_sum_factor
    ! type(cq_timer) :: tmr_l_tmp1, tmr_l_iter

    call start_timer(tmr_l_tmp1, WITH_LEVEL)
    prefix = return_prefix(subname, min_layer)
    spin_SF = 1
    if(inode==ionode .and. iprint_DM + min_layer >= 2) write(io_lun,fmt='(/4x,a)') trim(prefix)//" starting"
    ! allocate temp matrices
    do spin = 1, nspin
       matM3(spin) = allocate_temp_matrix(Lrange,0)
       matSM3(spin) = allocate_temp_matrix(Lrange,0)
       matSphi(spin) = allocate_temp_matrix(Lrange,0)
       mat_temp(spin) = allocate_temp_matrix(TLrange,0)
       mat_search(spin) = allocate_temp_matrix(Lrange,0)
    end do

    if (ndone > n_L_iterations) &
         call cq_abort('earlyDM: too many L iterations', &
                       ndone, n_L_iterations)

    min_layer = min_layer - 1
    if (vary_mu) call correct_electron_number
    min_layer = min_layer + 1

    call LNV_matrix_multiply(electrons, energy0, dontK, dontM1, &
                             dontM2, doM3, dontM4, dophi, doE,  &
                             mat_M3=matM3, mat_phi=matphi)
    electrons_tot = spin_factor * sum(electrons(:))
    energy0_tot = spin_factor * sum(energy0(:))

    do spin = 1, nspin
       if (flag_SpinDependentSF) spin_SF = spin
       ! Pre- and post-multiply M3 by S^-1 so that it is contravariant
       call matrix_product(matT(spin_SF), matM3(spin), mat_temp(spin), mult(T_L_TL))
       call matrix_product(mat_temp(spin), matT(spin_SF), matSM3(spin), mult(TL_T_L))
       ! Project gradient perpendicular to electron gradient
       call matrix_sum(zero, mat_search(spin), -one, matSM3(spin))
    end do

    if (vary_mu) then
       e_dot_n_tot = zero
       n_dot_n_tot = zero
       do spin = 1, nspin
          if (flag_SpinDependentSF) spin_SF = spin
          ! Pre- and post-multiply phi by S^-1 so that it is contravariant
          call matrix_product(matT(spin_SF), matphi(spin), mat_temp(spin), mult(T_L_TL))
          call matrix_product(mat_temp(spin), matT(spin_SF), matSphi(spin), mult(TL_T_L))
          ! Only one is pre- and post-multiplied because (A,B) = Tr(ASBS)
          e_dot_n(spin) = matrix_product_trace(matSM3(spin), matphi(spin))
          e_dot_n_tot = e_dot_n_tot + spin_factor * e_dot_n(spin)
          n_dot_n(spin) = matrix_product_trace(matSphi(spin), matphi(spin))
          n_dot_n_tot = n_dot_n_tot + spin_factor * n_dot_n(spin)
          if (inode == ionode .and. iprint_DM + min_layer >= 3) then
             write(io_lun, '(4x,a,i1,") ",f16.6)') &
                  trim(prefix)//" e.n (spin=", spin, e_dot_n(spin)
             write(io_lun, '(4x,a,i1,") ",f16.6)') &
                  trim(prefix)//" n.n (spin=", spin, n_dot_n(spin)
          end if
       end do
       ! Correct search direction so that it is tangent to
       ! iso-surface of electron numbers. This is different for if
       ! one fixes spin populations or not
       ! If spin population is fixed then (x denote spin)
       !
       !                     tr(dE/dL^x * dN/dL^x)
       !   G^x = - dE/dL^x + --------------------- dN/dL^x
       !                     tr(dN/dL^x * dN/dL^x)
       !
       ! If only the total electron number is fixed then
       !
       !                       sum_x tr(dE/dL^x * dN/dL^x)
       !   direct_sum_factor = ---------------------------
       !                       sum_x tr(dN/dL_x * dN/dL^x)
       !
       !   G^x = - dE/dL^x + direct_sum_factor * dN/dL^x
       !
       if (nspin == 1 .or. flag_fix_spin_population) then
          do spin = 1, nspin
             ! if gradient of N w.r.t. L is zero then no correction needed
             if (abs(n_dot_n(spin)) > RD_ERR) then
                call matrix_sum(one, mat_search(spin), &
                                e_dot_n(spin) / n_dot_n(spin), &
                                matSphi(spin))
             end if
          end do
       else ! variable spin
          ! if gradient of N w.r.t. L is zero then no correction needed
          if (abs(n_dot_n_tot) > RD_ERR) then
             direct_sum_factor = e_dot_n_tot / n_dot_n_tot
             do spin = 1, nspin
                call matrix_sum(one, mat_search(spin), direct_sum_factor,&
                                matSphi(spin))
             end do
          end if
       end if ! (nspin == 1 .or. flag_fix_spin_population)
    end if ! (vary_mu)

    call stop_print_timer(tmr_l_tmp1, "earlyDM - Preliminaries", &
                          IPRINT_TIME_THRES1)

    !-----------!
    ! MAIN LOOP !
    !-----------!
    ! delta_e is used for spin polarised calculations too
    delta_e = zero
    do n_iter = 1, n_L_iterations
       call start_timer(tmr_l_iter, WITH_LEVEL)
       if (flag_global_tolerance) then
          g0 = zero
          do spin = 1, nspin
             g0 = g0 + spin_factor * &
                  matrix_product_trace(matM3(spin), matSM3(spin))
          end do
       else
          g0 = zero
          do spin = 1, nspin
             g0 = g0 + spin_factor * &
                  matrix_product_trace(matM3(spin), matSM3(spin)) / &
                  real(ni_in_cell, double)
          end do
       end if
       if (inode == ionode .and. iprint_DM + min_layer >= 1) &
            write(io_lun, fmt='(4x,a,i3," Energy: ",f16.6," Ha Residual: ",e16.6)') &
            trim(prefix)//" iteration: ",n_iter, energy0_tot, g0

       ! minimise total E along search direction, updates energy1 =
       ! E(L_n_iter+1), delta_e = energy1_tot - energy0_tot,
       ! matM3(L_n_iter+1), matM3_dn(L_n_iter+1), inflex and interpG
       min_layer = min_layer - 1
       call lineMinL(matM3, mat_search, mat_temp, matSM3,   &
                     energy0_tot, energy1_tot, delta_e, &
                     inflex, interpG)
       min_layer = min_layer + 1
       ! panic if found inflexion point
       if (inflex) then
          call cq_warn(subname,"Inflexion point found in DM search; resetting")
          ndone = n_iter
          ! deallocate matrices
          do spin = nspin, 1, -1
             call free_temp_matrix(mat_search(spin))
             call free_temp_matrix(mat_temp(spin))
             call free_temp_matrix(matSphi(spin))
             call free_temp_matrix(matSM3(spin))
             call free_temp_matrix(matM3(spin))
          end do
          call stop_print_timer(tmr_l_iter, 'an earlyDM iteration', &
                                IPRINT_TIME_THRES1)
          return !This is a panic sign !
       end if ! inflex

       ! Gradient after line min - assumes LVN_matrix_mult at end
       ! of lineMinL
       ! Pre- and post-multiply M3 by S^-1 so that it is
       ! contravariant, matSM3 = matSM3(L_n_iter+1)
       do spin = 1, nspin
          if (flag_SpinDependentSF) spin_SF = spin
          call matrix_product (matT(spin_SF), matM3(spin), mat_temp(spin), mult(T_L_TL))
          call matrix_product (mat_temp(spin), matT(spin_SF), matSM3(spin), mult(TL_T_L))
       end do

       ! g1 = tr(matM3(L_n_iter+1) * matSM3(L_n_iter+1))
       if(flag_global_tolerance) then
          g1 = zero
          do spin = 1, nspin
             g1 = g1 + spin_factor * &
                  matrix_product_trace(matM3(spin), matSM3(spin))
          end do
       else
          g1 = zero
          do spin = 1, nspin
             g1 = g1 + spin_factor * &
                  matrix_product_trace(matM3(spin), matSM3(spin)) / &
                  real(ni_in_cell, double)
          end do
       end if

       ! Test for linearity or convergence
       ! zeta is returned by lineMinL
       zeta = interpG
       if (inode == ionode .and. iprint_DM + min_layer >= 2) &
            write(io_lun,fmt='(4x,a,f16.6,a,f16.6)') trim(prefix)//" dE ", delta_e, " Ha linearity: ",zeta
       if (inode == ionode .and. iprint_DM + min_layer >= 3) &
            write (io_lun, '(4x,a,3f16.6)') trim(prefix)//" zeta...: ", g0, g1, interpG

       ! Correct L_n_iter+1 again to get electron numbers correct
       min_layer = min_layer - 1
       if (vary_mu) call correct_electron_number
       min_layer = min_layer + 1

       ! recalculate the quantities after L_n_iter+1 is updated
       ! 2011/08/24 L.Tong:
       !   The orignal implementation uses energy0 as input for
       !   LNV_matrix_multiply, but I think this is incorrect and
       !   should be energy1
       call LNV_matrix_multiply(electrons, energy1, dontK, dontM1, &
                                dontM2, doM3, dontM4, dophi, doE, &
                                mat_M3=matM3, mat_phi=matphi)
       electrons_tot = spin_factor * sum(electrons(:))
       energy1_tot = spin_factor * sum(energy1(:))

       ! Pre- and post-multiply M3 by S^-1 so that it is contravariant
       do spin = 1, nspin
          if (flag_SpinDependentSF) spin_SF = spin
          call matrix_product(matT(spin_SF), matM3(spin), mat_temp(spin), mult(T_L_TL))
          call matrix_product(mat_temp(spin), matT(spin_SF), matSM3(spin), mult(TL_T_L))
       end do

       ! prepare for the next iterative step
       energy0_tot = energy1_tot
       ! update the search direction
       do spin = 1, nspin
          call matrix_sum(zero, mat_search(spin), -one, matSM3(spin))
       end do

       ! Project search direction perpendicular to electron gradient
       if (vary_mu) then
          ! Wrap electron gradient for tensorial correctness
          e_dot_n_tot = zero
          n_dot_n_tot = zero
          do spin = 1, nspin
             if (flag_SpinDependentSF) spin_SF = spin
             call matrix_product(matT(spin_SF), matphi(spin), mat_temp(spin), &
                                 mult(T_L_TL))
             call matrix_product(mat_temp(spin), matT(spin_SF), matSphi(spin), &
                                 mult(TL_T_L))
             e_dot_n(spin) = matrix_product_trace(matSM3(spin), matphi(spin))
             e_dot_n_tot = e_dot_n_tot + spin_factor * e_dot_n(spin)
             n_dot_n(spin) = matrix_product_trace(matSphi(spin), matphi(spin))
             n_dot_n_tot = n_dot_n_tot + spin_factor * n_dot_n(spin)
          end do
          if (nspin == 1 .or. flag_fix_spin_population) then
             do spin = 1, nspin
                ! if gradient of N w.r.t. L is zero then no correction needed
                if (abs(n_dot_n(spin)) > RD_ERR) then
                   call matrix_sum(one, mat_search(spin), &
                                   e_dot_n(spin) / n_dot_n(spin), &
                                   matSphi(spin))
                end if
             end do
          else
             ! if gradient of N w.r.t. L is zero then no correction needed
             if (abs(n_dot_n_tot) > RD_ERR) then
                direct_sum_factor = e_dot_n_tot / n_dot_n_tot
                do spin = 1, nspin
                   call matrix_sum(one, mat_search(spin), &
                                   direct_sum_factor, matSphi(spin))
                end do
             end if
          end if
          ! Here, we can't alter M3 directly: lineMinL expects REAL gradient
       end if

       if (flag_global_tolerance) then
          g0 = zero
          do spin = 1, nspin
             g0 = g0 + spin_factor * &
                  matrix_product_trace(matM3(spin), matSM3(spin))
          end do
       else
          g0 = zero
          do spin = 1, nspin
             g0 = g0 + spin_factor * &
                  matrix_product_trace(matM3(spin), matSM3(spin)) / &
                  real(ni_in_cell, double)
          end do
       end if

       ! store the Pulay histories
       PulayR(ndone + n_iter) = g0
       PulayE(ndone + n_iter) = energy0_tot

       ! test if we are linear
       if (abs(zeta) < LinTol_DMM) then
          ! We're linear Added +1 to n_iter for cosmetic reasons (this
          ! is true since at this stage g0 and energy0_tot etc are already
          ! prepared for n_iter+1 step)
          if (inode == ionode .and. iprint_DM + min_layer >= 1) &
               write(io_lun, fmt='(4x,a,i3," Energy: ",f16.6," Ha Residual: ",e16.6)') &
               trim(prefix)//" iteration: ",n_iter+1, energy0_tot, g0
          if ((inode == ionode).and. (iprint_DM + min_layer >= 1)) &
               write (io_lun, '(4x,a,i3,a)') &
               trim(prefix)//" transition to lateDM after ",n_iter+1," iterations"
          ndone = n_iter
          do spin = nspin, 1, -1
             call free_temp_matrix(mat_search(spin))
             call free_temp_matrix(mat_temp(spin))
             call free_temp_matrix(matSphi(spin))
             call free_temp_matrix(matSM3(spin))
             call free_temp_matrix(matM3(spin))
          end do
          call stop_print_timer(tmr_l_iter, "an earlyDM iteration", &
                                IPRINT_TIME_THRES1)
          return
       end if
       call stop_print_timer(tmr_l_iter, "an earlyDM iteration", &
                             IPRINT_TIME_THRES1)
    end do ! do n_iter = 1, n_L_iterations

    do spin = nspin, 1, -1
       call free_temp_matrix(mat_search(spin))
       call free_temp_matrix(mat_temp(spin))
       call free_temp_matrix(matSphi(spin))
       call free_temp_matrix(matSM3(spin))
       call free_temp_matrix(matM3(spin))
    end do
    return
  end subroutine earlyDM
  !!***


  ! -----------------------------------------------------------
  ! Subroutine lateDM
  ! -----------------------------------------------------------

  !!****f* DMMin/lateDM *
  !!
  !!  NAME
  !!   lateDM
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Performs late-stage (linear) minimisation of energy wrt density
  !!   matrix elements.  Uses the GR-Pulay scheme (described in Chem Phys Lett
  !!   325, 796 (2000)).  If the residual (given by the squared norm of the
  !!   gradient) increases, this is taken as a sign that linearity has been
  !!   broken, and the scheme exits to earlyDM
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   11/12/00
  !!  MODIFICATION HISTORY
  !!   05/06/2001 dave
  !!    Added ROBODoc header and removed stop commands (for cq_abort)
  !!   08/06/2001 dave
  !!    Added my_barrier to GenComms use list
  !!    Removed my_barrier call (unnecessary) and changed dgsum to gsum
  !!   2004/10/28 drb
  !!    Small tidying
  !!   10:09, 13/02/2006 drb
  !!    Removed all explicit references to data_ variables and rewrote
  !!    in terms of new matrix routines
  !!   2008/03/03 18:32 dave
  !!    Removed dsqrt
  !!   2009/07/08 16:41 dave
  !!    Introduced atom-based tolerance (just divide by ni_in_cell)
  !!   2010/03/18 14:08 dave
  !!    Added code for mixed L and SCF minimisation
  !!   2011/08/25 L.Tong
  !!    - Added spin polarisation
  !!    - NOTE:
  !!        There is a change to the size of residual. The new
  !!        residuals are now a factor 2 less than the original
  !!        implementation. The actual method for calculating the
  !!        residual has not changed. However since the definition of
  !!        mat_M3 has changed, now mat_M3(nspin) is an array and
  !!        mat_M3(1) always stores the values for spin up component,
  !!        this has resulted a changed in the calculated residual.
  !!
  !!        For spin non-polarised calculations:
  !!        The original implementation has mat_M3 = 2.0 * mat_M3(1)
  !!        Hence:  g1 = tr(mat_M3*T*mat_M3*T)
  !!                   = 4.0 * tr(mat_M3(1)*T*mat_M3(1)*T)
  !!        Now this has becomes
  !!                g1 = sum_{is} tr(mat_M3(is)*T*mat_M3(is)*T)
  !!                   = 2.0 * tr(mat_M3(1)*T*mat_M3(1)*T)
  !!        Hence the factor two difference.
  !!    - Removed redundant parameter number_of_bands
  !!   2012/03/12 L.Tong
  !!    - Major rewrite for change in spin implementation
  !!    - Removed redundant input paramter mu (real, double)
  !!   2012/06/24 L.Tong
  !!    - Added control of whether to dump L using flag_dump_L
  !!   2013/03/06 L.Tong
  !!    - Changed the min|max step size to minpulaystepDMM and maxpulaystepDMM
  !!      instead of hardwired values of 0.001 and 0.1
  !!   2013/08/20 M.Arita
  !!    - Changed calls for dumping L-matrix
  !!   2016/07/13 18:30 nakata
  !!    Renamed H_on_supportfns -> H_on_atomfns
  !!   2016/08/08 15:30 nakata
  !!    Renamed supportfns -> atomfns
  !!   2017/05/09 dave
  !!    Added output of L matrix for both spin channels
  !!   2018/11/13 17:30 nakata
  !!    Changed matT to be spin_SF dependent
  !!  SOURCE
  !!
  subroutine lateDM(ndone, n_L_iterations, done, deltaE, vary_mu, &
                    tolerance, record)

    use datatypes
    use numbers
    use logicals
    use matrix_data,       only: Lrange, TLrange
    use mult_module,       only: matrix_product, allocate_temp_matrix, &
                                 free_temp_matrix,                     &
                                 LNV_matrix_multiply, matphi, matT,    &
                                 matL, mult, T_L_TL, TL_T_L,           &
                                 matrix_sum, matrix_product_trace,     &
                                 matrix_scale
    use Pulay
    use PosTan,            only: PulayR, PulayE, max_iters
    use GenComms,          only: cq_abort, gsum
    use global_module,     only: IPRINT_TIME_THRES1,        &
                                 ni_in_cell, flag_global_tolerance,    &
                                 flag_mix_L_SC_min,                    &
                                 flag_fix_spin_population, nspin,      &
                                 spin_factor, flag_dump_L,             &
                                 flag_SpinDependentSF, min_layer
    use timer_module,      only: cq_timer,start_timer,                 &
                                 stop_print_timer, WITH_LEVEL
    use io_module,         only: dump_matrix, return_prefix
    use functions_on_grid, only: atomfns, H_on_atomfns
    use H_matrix_module,   only: get_H_matrix
    use density_module,    only: density, get_electronic_density, flag_DumpChargeDensity
    use maxima_module,     only: maxngrid
    use memory_module,     only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    !Prints out charge density -- 2010.Nov.06 TM
    use io_module,         only: dump_charge
    use dimens,            only: n_my_grid_points
    use store_matrix,      only: dump_pos_and_matrices, unit_DM_save

    implicit none

    ! Passed variables
    integer      :: n_L_iterations, length, ndone
    logical      :: vary_mu, done, record
    real(double) :: tolerance, deltaE

    ! Local variables
    integer, dimension(maxpulayDMM,nspin) :: mat_Lstore, mat_Gstore, mat_SGstore
    integer, dimension(nspin) :: matM3, matSM3, matSphi, mat_temp
    real(double), dimension(nspin) :: e_dot_n, n_dot_n, electrons, &
                                      energy0, energy1
    real(double), dimension(maxpulayDMM)             :: alph
    real(double), dimension(maxpulayDMM,maxpulayDMM) :: Aij
    ! real(double), dimension(maxpulayDMM*maxpulayDMM) :: Aij1
    real(double)   :: e_dot_n_tot, n_dot_n_tot, electrons_tot, &
                      energy0_tot, energy1_tot
    real(double)   :: g0, g1, gg, step, dsum_factor
    integer        :: n_iter, i, j, pul_mx, npmod, spin, spin_SF, stat
    type(cq_timer) :: tmr_l_tmp1,tmr_l_iter
    real(double), dimension(:), allocatable :: density_tot
    !TM
    integer        :: iter_stuck = 0
    integer, parameter :: mx_stuck = 5
    character(len=12) :: subname = "lateDM: "
    character(len=120) :: prefix

    call start_timer(tmr_l_tmp1, WITH_LEVEL)
    prefix = return_prefix(subname, min_layer)
    iter_stuck = 0

    spin_SF = 1

    if (ndone > n_L_iterations) &
         call cq_abort('lateDM: too many L iterations', ndone, n_L_iterations)

    do spin = 1, nspin
       do i = 1, maxpulayDMM
          mat_Lstore(i,spin) = allocate_temp_matrix(Lrange,0)
          mat_Gstore(i,spin) = allocate_temp_matrix(Lrange,0)
          mat_SGstore(i,spin) = allocate_temp_matrix(Lrange,0)
       end do
       matM3(spin) = allocate_temp_matrix(Lrange,0)
       matSM3(spin) = allocate_temp_matrix(Lrange,0)
       matSphi(spin) = allocate_temp_matrix(Lrange,0)
       mat_temp(spin) = allocate_temp_matrix(TLrange,0)
    end do

    ! Update the charge density if flag is set
    min_layer = min_layer - 1
    if (flag_mix_L_SC_min) then
       ! 2011/08/29 L.Tong:
       ! original also calculates matphi, (dophi), but I think
       ! this is redundant. So set dontphi. Only doK.
       call LNV_matrix_multiply(electrons, energy1, doK, dontM1, &
                                dontM2, dontM3, dontM4, dontphi, doE)
       energy1_tot = spin_factor * sum(energy1(:))

       ! note H_on_atomfns(1) is used just as a temp working array
       call get_electronic_density(density, electrons, atomfns,    &
                                   H_on_atomfns(1), inode, ionode, &
                                   maxngrid)
       call get_H_matrix(.true., .false., electrons, density, maxngrid)
    end if
    min_layer = min_layer + 1

    ! Get the gradient at the starting point (?) this updates matM3
    ! and matphi
    call LNV_matrix_multiply(electrons, energy0, dontK, dontM1, &
                             dontM2, doM3, dontM4, dophi, doE,  &
                             mat_M3=matM3, mat_phi=matphi)
    electrons_tot = spin_factor * sum(electrons)
    energy0_tot = spin_factor * sum(energy0)

    ! Covariant gradient in SM3
    do spin = 1, nspin
       if (flag_SpinDependentSF) spin_SF = spin
       call matrix_product(matT(spin_SF), matM3(spin), mat_temp(spin), mult(T_L_TL))
       call matrix_product(mat_temp(spin), matT(spin_SF), matSM3(spin), mult(TL_T_L))
    end do

    ! Project electron gradient out
    if (vary_mu) then
       ! update matSphi from matphi
       do spin = 1, nspin
          if (flag_SpinDependentSF) spin_SF = spin
          call matrix_product(matT(spin_SF), matphi(spin), mat_temp(spin), mult(T_L_TL))
          call matrix_product(mat_temp(spin), matT(spin_SF), matSphi(spin), mult(TL_T_L))
          e_dot_n(spin) = matrix_product_trace(matSM3(spin), matphi(spin))
          n_dot_n(spin) = matrix_product_trace(matSphi(spin), matphi(spin))
          if (inode == ionode .and. iprint_DM + min_layer >= 3) then
             write(io_lun, '(4x,a,i1,") ",f16.6)') &
                  trim(prefix)//" e.n (spin=", spin, e_dot_n(spin)
             write(io_lun, '(4x,a,i1,") ",f16.6)') &
                  trim(prefix)//" n.n (spin=", spin, n_dot_n(spin)
          end if
       end do

       ! matSM3 mat matM3 are used as search direction (-SG and -G)
       ! These should be positive because we SUBTRACT M3 off below
       ! Here, we can alter M3 directly because it's not expected
       ! to be an exact gradient
       if (nspin == 1 .or. flag_fix_spin_population) then
          do spin = 1, nspin
             ! if gradient of N w.r.t. L is zero then no correction needed
             if (abs(n_dot_n(spin)) > RD_ERR) then
                dsum_factor = e_dot_n(spin) / n_dot_n(spin)
                call matrix_sum(one, matSM3(spin), -dsum_factor, matSphi(spin))
                call matrix_sum(one, matM3(spin), -dsum_factor, matphi(spin))
             end if
          end do
       else
          e_dot_n_tot = spin_factor * sum(e_dot_n)
          n_dot_n_tot = spin_factor * sum(n_dot_n)
          ! if gradient of N w.r.t. L is zero then no correction needed
          if (abs(n_dot_n_tot) > RD_ERR) then
             dsum_factor = e_dot_n_tot / n_dot_n_tot
             do spin = 1, nspin
                call matrix_sum(one, matSM3(spin), -dsum_factor, matSphi(spin))
                call matrix_sum(one, matM3(spin), -dsum_factor, matphi(spin))
             end do
          end if
       end if
    end if ! (vary_mu)
    ! Store initial gradient (these are search directions G and L)
    do spin = 1, nspin
       call matrix_sum(zero, mat_SGstore(1,spin), -one, matSM3(spin))
       call matrix_sum(zero, mat_Gstore(1,spin), -one, matM3(spin))
       call matrix_sum (zero, mat_Lstore(1,spin), one, matL(spin))
    end do

    ! timer
    call stop_print_timer (tmr_l_tmp1, "lateDM - Preliminaries", &
                           IPRINT_TIME_THRES1)

    !-----------!
    ! MAIN LOOP !
    !-----------!

    if (flag_global_tolerance) then
       g0 = zero
       do spin = 1, nspin
          g0 = g0 + spin_factor * &
               matrix_product_trace(matM3(spin), matSM3(spin))
       end do
    else
       g0 = zero
       do spin = 1, nspin
          g0 = g0 + spin_factor * &
               matrix_product_trace(matM3(spin), matSM3(spin)) / &
               real(ni_in_cell, double)
       end do
    end if
    if(n_L_iterations==0) done = .true. ! Otherwise the calling routine will cycle
    do n_iter = 1, n_L_iterations
       call start_timer(tmr_l_iter, WITH_LEVEL)
       if (inode == ionode .and. iprint_DM + min_layer >= 1) &
            write(io_lun, fmt='(4x,a,i3," Energy: ",f16.6," Ha Residual: ",e16.6)') &
            trim(prefix)//"  iteration: ",n_iter, energy0_tot, g0
       ! Storage for pulay DMs/residuals
       npmod = mod(n_iter, maxpulayDMM) + 1
       pul_mx = min(n_iter + 1, maxpulayDMM)
       ! Take a step - maybe correct electron number after
       if (flag_global_tolerance) then
          ! Base step on present gradient and expected dE
          step = deltaE / (g0 + RD_ERR)
       else
          ! Base step on present gradient and expected dE
          step = deltaE / (real(ni_in_cell, double) * g0 + RD_ERR)
       end if
       step = abs(step)
       ! We don't want the step to be too small or too big
       if (abs(step) < minpulaystepDMM) step = minpulaystepDMM
       if (abs(step) > maxpulaystepDMM) step = maxpulaystepDMM
       if (inode == ionode .and. iprint_DM + min_layer >= 3) &
            write(io_lun, '(4x,a,i3,i3,f16.6)') &
            trim(prefix)//" npmod, pul_mx and step: ", npmod, pul_mx, step
       ! take L_n+1 = L_n + step * G_n
       if (npmod > 1) then
          do spin = 1, nspin
             call matrix_sum(zero, matL(spin), one, mat_Lstore(npmod-1,spin))
             call matrix_sum(one, matL(spin), step, mat_SGstore(npmod-1,spin))
          end do
       else
          do spin = 1, nspin
             call matrix_sum(zero, matL(spin), one, mat_Lstore(pul_mx,spin))
             call matrix_sum(one, matL(spin), step, mat_SGstore(pul_mx,spin))
          end do
       endif
       ! after the step, correct the electron number
       min_layer = min_layer - 1
       if (vary_mu) call correct_electron_number
       min_layer = min_layer + 1
       ! Re-evaluate gradient and energy
       call LNV_matrix_multiply(electrons, energy1, dontK, dontM1, &
                                dontM2, doM3, dontM4, dophi, doE,  &
                                mat_M3=matM3, mat_phi=matphi)
       electrons_tot = spin_factor * sum(electrons)
       energy1_tot = spin_factor * sum(energy1)
       ! Covariant gradient in SM3
       do spin = 1, nspin
          if (flag_SpinDependentSF) spin_SF = spin
          call matrix_product(matT(spin_SF), matM3(spin), mat_temp(spin), mult(T_L_TL))
          call matrix_product(mat_temp(spin), matT(spin_SF), matSM3(spin), mult(TL_T_L))
       end do
       ! Project out electron variation
       if (vary_mu) then
          ! update matSphi
          do spin = 1, nspin
             if (flag_SpinDependentSF) spin_SF = spin
             call matrix_product(matT(spin_SF), matphi(spin), mat_temp(spin), &
                                 mult(T_L_TL))
             call matrix_product(mat_temp(spin), matT(spin_SF), matSphi(spin), &
                                 mult(TL_T_L))
             e_dot_n(spin) = matrix_product_trace(matSM3(spin), matphi(spin))
             n_dot_n(spin) = matrix_product_trace(matSphi(spin), matphi(spin))
          end do
          if (flag_fix_spin_population) then
             do spin = 1, nspin
                ! if gradient of N w.r.t. L is zero then no correction needed
                if (abs(n_dot_n(spin)) > RD_ERR) then
                   dsum_factor = e_dot_n(spin) / n_dot_n(spin)
                   call matrix_sum(one, matSM3(spin), -dsum_factor, &
                                   matSphi(spin))
                   call matrix_sum(one, matM3(spin), -dsum_factor, &
                                   matphi(spin))
                end if
             end do
          else
             e_dot_n_tot = spin_factor * sum(e_dot_n)
             n_dot_n_tot = spin_factor * sum(n_dot_n)
             ! if gradient of N w.r.t. L is zero then no correction needed
             if (abs(n_dot_n_tot) > RD_ERR) then
                dsum_factor = e_dot_n_tot / n_dot_n_tot
                do spin = 1, nspin
                   call matrix_sum (one, matSM3(spin), -dsum_factor, &
                                    matSphi(spin))
                   call matrix_sum (one, matM3(spin), -dsum_factor, &
                                    matphi(spin))
                end do
             end if
          end if
       end if
       ! 2022/10/28 16:01 dave
       ! I'm not sure that we need this calculation of residual
       ! Find the residual (i.e. the gradient)
       if (flag_global_tolerance) then
          gg = zero
          do spin = 1, nspin
             gg = gg + spin_factor * &
                  matrix_product_trace(matM3(spin), matSM3(spin))
          end do
       else
          gg = zero
          do spin = 1, nspin
             gg = gg + spin_factor * &
                  matrix_product_trace(matM3(spin), matSM3(spin)) / &
                  real(ni_in_cell, double)
          end do
       end if
       if (inode == ionode .and. iprint_DM + min_layer >= 2) &
            write(io_lun, '(4x,a,e16.6)') trim(prefix)//" R2 is ", sqrt(gg)
       ! record Pulay histories
       do spin = 1, nspin
          call matrix_sum(zero, mat_SGstore(npmod,spin), -one, matSM3(spin))
          call matrix_sum(zero, mat_Gstore(npmod,spin), -one, matM3(spin))
          call matrix_sum(zero, mat_Lstore(npmod,spin), one, matL(spin))
       end do
       ! calculate A_ij
       Aij = zero
       do i = 1, pul_mx
          do j = 1, pul_mx
             gg = zero
             do spin = 1, nspin
                gg = gg + spin_factor * &
                     matrix_product_trace(mat_Gstore(j,spin),&
                                          mat_SGstore(i,spin))
             end do
             Aij(j,i) = gg
             ! Aij1(j + (i-1)*pul_mx) = gg
          end do ! j
       end do ! i
       ! Solve to get alphas
       call DoPulay(npmod, Aij, alph, pul_mx, maxpulayDMM)
       ! Make new L matrix from Pulay sum
       do spin = 1, nspin
          call matrix_scale(zero, matL(spin))
          do i = 1, pul_mx
             call matrix_sum (one, matL(spin), alph(i), mat_Lstore(i,spin))
          end do
       end do
       ! after the step, correct the electron number
       min_layer = min_layer - 1
       if (vary_mu) call correct_electron_number
       if (flag_mix_L_SC_min) then
          ! 2011/08/29 L.Tong
          ! the original also has dophi, but I think this is redundant
          call LNV_matrix_multiply (electrons, energy1, doK, dontM1, &
                                    dontM2, dontM3, dontM4, dontphi, doE)
          ! H_on_atomfns(1) is used just as a temp working array
          call get_electronic_density(density, electrons,     &
                                      atomfns,                &
                                      H_on_atomfns(1), inode, &
                                      ionode, maxngrid)
          call get_H_matrix(.true., .false., electrons, density, maxngrid)
       end if
       min_layer = min_layer + 1
       ! re-evaluate the gradient and energy at new position
       call LNV_matrix_multiply(electrons, energy1, dontK, dontM1, &
                                dontM2, doM3, dontM4, dophi, doE,  &
                                mat_M3=matM3, mat_phi=matphi)
       energy1_tot = spin_factor * sum(energy1)
       electrons_tot = spin_factor * sum(electrons)
       do spin = 1, nspin
          if (flag_SpinDependentSF) spin_SF = spin
          call matrix_product(matT(spin_SF), matM3(spin), mat_temp(spin), mult(T_L_TL))
          call matrix_product(mat_temp(spin), matT(spin_SF), matSM3(spin), mult(TL_T_L))
       end do
       if (flag_global_tolerance) then
          g1 = zero
          do spin = 1, nspin
             g1 = g1 + spin_factor * &
                  matrix_product_trace(matM3(spin), matSM3(spin))
          end do
       else
          g1 = zero
          do spin = 1, nspin
             g1 = g1 + spin_factor * &
                  matrix_product_trace(matM3(spin), matSM3(spin)) / &
                  real(ni_in_cell, double)
          end do
       end if
       if (inode == ionode .and. iprint_DM + min_layer >= 2) &
            write(io_lun, '(4x,a,e16.6)') &
            trim(prefix)//" Residual before electron gradient correction: ", g1
       if (vary_mu) then
          do spin = 1, nspin
             if (flag_SpinDependentSF) spin_SF = spin
             call matrix_product(matT(spin_SF), matphi(spin), &
                                 mat_temp(spin), mult(T_L_TL))
             call matrix_product(mat_temp(spin), matT(spin_SF), &
                                 matSphi(spin), mult(TL_T_L))
             e_dot_n(spin) = matrix_product_trace(matSM3(spin), matphi(spin))
             n_dot_n(spin) = matrix_product_trace(matSphi(spin), matphi(spin))
          end do
          if (flag_fix_spin_population) then
             do spin = 1, nspin
                ! if gradient of N w.r.t. L is zero then no correction needed
                if (abs(n_dot_n(spin)) > RD_ERR) then
                   dsum_factor = e_dot_n(spin) / n_dot_n(spin)
                   call matrix_sum(one, matSM3(spin), -dsum_factor, &
                                   matSphi(spin))
                   call matrix_sum(one, matM3(spin), -dsum_factor, &
                                   matphi(spin))
                end if
             end do
          else
             e_dot_n_tot = spin_factor * sum(e_dot_n)
             n_dot_n_tot = spin_factor * sum(n_dot_n)
             ! if gradient of N w.r.t. L is zero then no correction needed
             if (abs(n_dot_n_tot) > RD_ERR) then
                dsum_factor = e_dot_n_tot / n_dot_n_tot
                do spin = 1, nspin
                   call matrix_sum(one, matSM3(spin), -dsum_factor, &
                                   matSphi(spin))
                   call matrix_sum(one, matM3(spin), -dsum_factor, &
                                   matphi(spin))
                end do
             end if
          end if
       end if
       ! get deltaE
       deltaE = energy1_tot - energy0_tot
       dE_DMM = deltaE
       ! Find the residual
       if (flag_global_tolerance) then
          g1 = zero
          do spin = 1, nspin
             g1 = g1 + spin_factor * &
                  matrix_product_trace(matM3(spin), matSM3(spin))
          end do
       else
          g1 = zero
          do spin = 1, nspin
             g1 = g1 + spin_factor * &
                  matrix_product_trace(matM3(spin), matSM3(spin)) / &
                  real(ni_in_cell, double)
          end do
       end if
       if (inode == ionode .and. iprint_DM + min_layer >= 2) &
            write(io_lun, '(4x,a,e16.6)') trim(prefix)//" New residual: ", g1
       if ((ndone + n_iter) < max_iters) then
          PulayR(ndone + n_iter) = g1
          PulayE(ndone + n_iter) = energy1_tot
       end if

       ! dump the L matrix if required
       if (n_dumpL>0) then
          if(mod (n_iter, n_dumpL) == 0) call dump_pos_and_matrices(index=unit_DM_save)
       endif
       !if (flag_dump_L) then
       !   if (mod (n_iter, n_dumpL) == 0) then
       !    call dump_pos_and_matrices
       !   end if
       !end if

       ! check if tolerance is reached
       if (g1 < tolerance) then
          done = .true.
          ndone = n_iter
          if (inode == ionode .and. iprint_DM + min_layer >= 0) &
               write(io_lun, '(/4x,a,e16.6,a,i3,a/)') &
               trim(prefix)//" reached residual of ", g1, " after ",n_iter, " iterations"
          call dealloc_lateDM
          call stop_print_timer(tmr_l_iter, "a lateDM iteration", &
                                IPRINT_TIME_THRES1)
          return
       else if ((.not. flag_mix_L_SC_min) .and. (g1 > two * g0)) then
          if (inode == ionode .and. iprint_DM + min_layer >= 0) &
               write(io_lun, '(/4x,a,i3,a,f16.6,a,f16.6)') &
               trim(prefix)//" residual rise after  ", n_iter, &
               " iterations with energy and residual: ", energy1_tot, " Ha ", g1
          ndone = n_iter
          call dealloc_lateDM
          call stop_print_timer(tmr_l_iter, "a lateDM iteration", &
               IPRINT_TIME_THRES1)
          return
       else if (g1 > 0.99_double * g0) then
          iter_stuck = iter_stuck + 1
          if (iter_stuck > mx_stuck) then
             done = .false.
             ndone = n_iter
             if (inode == ionode .and. iprint_DM + min_layer >= 0) &
                  write(io_lun, '(/4x,a,i3,a,f16.6,a,e16.6)') &
                  trim(prefix)//" reset since failed to reduce residual after ",n_iter, &
                  " iterations with energy and residual: ", energy1_tot, " Ha ", g1
             call dealloc_lateDM
             call stop_print_timer(tmr_l_iter, "a lateDM iteration", &
                                   IPRINT_TIME_THRES1)
             return
          end if ! (iter_stuck > mx_stuck)
       end if
       ! Replace step with that calculated from real L, and prepare
       ! for next step
       do spin = 1, nspin
          call matrix_sum(zero, mat_SGstore(npmod,spin), -one, matSM3(spin))
          call matrix_sum(zero, mat_Gstore(npmod,spin), -one, matM3(spin))
          call matrix_sum(zero, mat_Lstore(npmod,spin), one, matL(spin))
       end do
       g0 = g1
       energy0_tot = energy1_tot
       call stop_print_timer(tmr_l_iter, "a lateDM iteration", &
                             IPRINT_TIME_THRES1)
    end do

    !Prints out charge density when selected
    if (flag_DumpChargeDensity .and. flag_mix_L_SC_min) then
       if (nspin == 1) then
          allocate(density_tot(maxngrid), STAT=stat)
          if (stat /= 0) call cq_abort("Error allocating density_tot: ", maxngrid)
          call reg_alloc_mem(area_DM, maxngrid, type_dbl)
          density_tot = zero
          density_tot = spin_factor * density(:,1)
          call dump_charge(density_tot, n_my_grid_points, inode)
          deallocate(density_tot, STAT=stat)
          if (stat /= 0) call cq_abort("Error deallocating density_tot")
          call reg_dealloc_mem(area_DM, maxngrid, type_dbl)

       else
          call dump_charge(density(:,1), n_my_grid_points, inode, spin=1)
          call dump_charge(density(:,2), n_my_grid_points, inode, spin=2)
       end if
    endif  ! (flag_mix_L_SC_min) then

    call dealloc_lateDM

    return

  contains

    subroutine dealloc_lateDM
      do spin = nspin, 1, -1
         call free_temp_matrix(mat_temp(spin))
         call free_temp_matrix(matSphi(spin))
         call free_temp_matrix(matSM3(spin))
         call free_temp_matrix(matM3(spin))
         do i = maxpulayDMM, 1, -1
            call free_temp_matrix(mat_SGstore(i,spin))
            call free_temp_matrix(mat_Gstore(i,spin))
            call free_temp_matrix(mat_Lstore(i,spin))
         end do
      end do
      return
    end subroutine dealloc_lateDM

  end subroutine lateDM
  !!***



! -----------------------------------------------------------
! Subroutine lineMinL
! -----------------------------------------------------------

  !!****f* DMMin/lineMinL_nospin *
  !!
  !!  NAME
  !!   lineMinL_nospin
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Performs an analytical line minimisation in the direction given by
  !!   mat_D.  Uses energy and gradient at 0 and step (where step is given by
  !!   the previous energy change) as the four pieces of information needed to
  !!   get the cubic coefficients analytically.  If the cubic has imaginary
  !!   maxima and minima, the solution is a point of inflexion, and we panic
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   D.R.Bowler (based on original by CMG and EHH)
  !!  CREATION DATE
  !!   11/12/00
  !!  MODIFICATION HISTORY
  !!   05/06/2001 dave
  !!    Added ROBODoc header
  !!   08/06/2001 dave
  !!    Changed dgsum to gsum
  !!   2004/10/28 drb
  !!    Fixed (partly ?) definition of linearity and changed to return
  !!    zeta in interpG
  !!   10:09, 13/02/2006 drb
  !!    Removed all explicit references to data_ variables and
  !!    rewrote in terms of new matrix routines
  !!   2011/08/23 L.Tong
  !!    Fixed a possible bug, if A=0, then the solution for true step should be
  !!    truestep = -C/2B instead of C/2B. Because we are solving the equation
  !!    3Ax^2 + 2Bx + C = 0, and if A = 0, x = -C/2B.
  !!   2012/03/01 L.Tong
  !!    renamed to lineMinL_nospin
  !!   2012/03/12 L.Tong
  !!   - Major rewrite to include spin implmentation. Now spin
  !!     polarised and non-polarised calculations share the same
  !!     subroutine.
  !!   - Renamed to lineMinL
  !!   2018/11/13 17:30 nakata
  !!    Changed matT to be spin_SF dependent
  !!  SOURCE
  !!
  subroutine lineMinL( matM3, mat_D, mat_temp, matSM3,  &
                      energy_in, energy_out, delta_e, &
                      inflex, interpG)

    use datatypes
    use numbers
    use logicals
    use matrix_data,   only: Lrange
    use mult_module,   only: matrix_sum, matrix_product_trace,       &
                             allocate_temp_matrix, free_temp_matrix, &
                             matrix_product, matT, matL,             &
                             LNV_matrix_multiply, mult, T_L_TL,      &
                             TL_T_L, symmetrise_L
    use global_module, only: flag_fix_spin_population, nspin,        &
                             spin_factor, flag_SpinDependentSF, min_layer
    use io_module,      only: return_prefix

    implicit none

    ! Passed variables
    integer, dimension(:) :: matM3, matSM3, mat_D, mat_temp
    logical               :: inflex
    real(double)          :: delta_e, energy_in, energy_out, interpG

    ! Local variables
    integer,      dimension(nspin) :: matM3old
    real(double), dimension(nspin) :: electrons_spin, energy_1_spin, energy_2_spin
    real(double) :: electrons, energy_0, energy_1
    real(double) :: step, truestep, directsum_step
    real(double) :: A, B, C, D, SQ
    real(double) :: g0, g1, ig0, ig1, ig1_step, igcross, igcross2
    real(double) :: lamG, zeta
    integer      :: i, j, k, length, spin, spin_SF
    character(len=12) :: subname = "lineMinL: "
    character(len=120) :: prefix

    prefix = return_prefix(subname, min_layer)
    spin_SF = 1
    g0 = zero
    ig0 = zero
    do spin = 1, nspin
       matM3old(spin) = allocate_temp_matrix(Lrange, 0)
       ! copy matM3(L_old)
       call matrix_sum(zero, matM3old(spin), one, matM3(spin))
       ! calculate g0 = tr(matD * matM3), trace of direct sum is sum of
       ! traces of spin components
       g0 = g0 + spin_factor * &
            matrix_product_trace(mat_D(spin), matM3(spin))

       ! calculate ig0
       ! L.Tong: WARNING! I just copied the approach from
       ! lineMinL_nospin, but I think this may be a mistake, and should
       ! be: ig0 = tr(matSM3 * matM3), not tr(mat_D, matM3)

       ig0 = ig0 + spin_factor * &
             matrix_product_trace(mat_D(spin), matM3(spin))
       ! ig0 = ig0 + spin_factor * &
       !       matrix_product_trace(matSM3(spin), matM3(spin))
    end do

    ! We have a cubic in energy with L_out = L + x.D
    ! energy at x = 0 is energy_in
    energy_0 = energy_in

    ! The intermediate step size should depend on expected energy reduction
    step = delta_e / g0
    if (abs(step) < 1.0e-2_double) step = 0.01_double
    if (abs(step) > 0.1_double)    step = 0.1_double
    if (inode == ionode .and. iprint_DM + min_layer >= 2) &
         write(io_lun, '(4x,a,f16.6)') trim(prefix)//" trial step is ",step

    ! now we take the step in the direction D
    do spin = 1, nspin
       call matrix_sum(one, matL(spin), step, mat_D(spin))
    end do

    ! calcuate matM3, energy and numbers of electrons again for L_step
    call LNV_matrix_multiply(electrons_spin, energy_1_spin, doK,    &
                             dontM1, dontM2, doM3, dontM4, dontphi, &
                             doE, mat_M3=matM3)
    energy_1 = spin_factor * sum(energy_1_spin(:))

    ! g1 = tr(mat_D(L_old) * matM3(L_step))
    g1 = zero
    do spin = 1, nspin
       g1 = g1 + spin_factor * &
            matrix_product_trace(mat_D(spin), matM3(spin))
    end do

    ! Evaluate old and new gradient cross and new gradient
    ! magnitude. This is for linear interpolation.
    ! igcross = tr(matSM3(L_old) * matM3(L_step))
    igcross = zero
    do spin = 1, nspin
       igcross = igcross + spin_factor * &
                 matrix_product_trace(matSM3(spin), matM3(spin))
    end do

    ! get the matSM3(L_step)
    do spin = 1, nspin
       if (flag_SpinDependentSF) spin_SF = spin
       call matrix_product(matT(spin_SF), matM3(spin), mat_temp(spin), mult(T_L_TL))
       call matrix_product(mat_temp(spin), matT(spin_SF), matSM3(spin), mult(TL_T_L))
    end do

    ! get ig1 = tr(matSM3(L_step) * matM3(L_step)),
    ! this is for linear interpolation
    ig1 = zero
    do spin = 1, nspin
       ig1 = ig1 + spin_factor * &
             matrix_product_trace(matSM3(spin), matM3(spin))
    end do
    ! remember the content of ig1 as it will be used to store other values
    ig1_step = ig1

    ! the cubic is then given by Ax^3+Bx^2+Cx+D with
    D = energy_0
    C = g0
    B = three * (energy_1 - D) / (step * step) - (g1 - two * C) / step
    A = (g1 - two * B * step - C) / (three * step * step)
    ! check for -ve square root
    SQ = four * B * B - 12.0_double * A * C

    if (SQ < 0) then
       ! point of inflexion - problem !
       if (inode == ionode .and. iprint_DM + min_layer > 0) then
          write(io_lun, fmt='(4x,a)') trim(prefix)//' inflexion approximation:'
          write(io_lun, fmt='(4x,a,4f16.6)') trim(prefix)//' A, B, C, D: ', A, B, C, D
       end if
       inflex = .true.
       do spin = nspin, 1, -1
          call free_temp_matrix(matM3old(spin))
       end do
       return
    else ! Solve cubic
       if (abs(A) > RD_ERR) then
          truestep = (-two * B + SQRT(SQ)) / (six * A)
       else
          if (inode == ionode .and. iprint_DM + min_layer >2) &
               write(io_lun, fmt='(4x,a)') trim(prefix)//' local quadratic approximation'
          ! L.Tong shouldn this be -C/(two*B) instead of the orginal
          ! C/(two * B) ?
          ! truestep = C / (two * B)
          truestep = - C / (two * B)
       end if
    end if ! if (SQ < 0) then
    if (inode == ionode .and. iprint_DM + min_layer>= 2) &
         write(io_lun, '(4x,a,f16.6)') trim(prefix)//" actual step is ", truestep

    !TM 09/09/2003
    if (truestep < 1.e-04_double) truestep = 1.e-04_double
    !TM 09/09/2003

    ! update value of matL
    directsum_step = truestep - step
    do spin = 1, nspin
       call matrix_sum(one, matL(spin), directsum_step, mat_D(spin))
    end do

    ! matrix symmetric: this avoids the creeping in of asymmetries due
    ! to accumulation of errors
    call symmetrise_L()

    ! Get K(L_truestep), E(L_truestep) and M3(L_truestep) before we return
    call LNV_matrix_multiply(electrons_spin, energy_2_spin, doK,    &
                             dontM1, dontM2, doM3, dontM4, dontphi, &
                             doE, mat_M3=matM3)
    energy_out = spin_factor * sum(energy_2_spin(:))
    delta_e = energy_out - energy_in

    !-----------------------------------------------------------------------
    !  finished calculating energy_out, delta_e, matM3, inflex
    !-----------------------------------------------------------------------

    ! Find interpolated gradient for linearity test (find interpG/zeta)
    if (truestep < step) then
       ! We need to compare interpG with g at truestep
       lamG = truestep / step
       interpG = ig0 * (one - lamG) * (one - lamG) + ig1 * lamG * lamG + &
                 two * igcross * lamG * (one - lamG)
       zeta = (interpG - ig1) / (ig0 - ig1)
    else ! truestep is GREATER than step
       ! We need to compare interpG with g at step
       ! find matSM3(L_truestep)
       ig1 = zero
       igcross = zero
       do spin = 1, nspin
          if (flag_SpinDependentSF) spin_SF = spin
          call matrix_product(matT(spin_SF), matM3(spin), mat_temp(spin), mult(T_L_TL))
          call matrix_product(mat_temp(spin), matT(spin_SF), matSM3(spin), mult(TL_T_L))
          ! ig1 is now used to store tr(matSM3(L_truestep) * matM3(L_truestep))
          ! ig1_step stores tr(matSM3(L_step) * matM3(L_step))
          ig1 = ig1 + spin_factor * &
                matrix_product_trace(matSM3(spin), matM3(spin))
          igcross = igcross + spin_factor * &
                    matrix_product_trace(matSM3(spin), matM3old(spin))
       end do
       lamG = step / truestep
       interpG = ig0 * (one - lamG) * (one - lamG) + ig1 * lamG * lamG + &
                 two * igcross * lamG * (one - lamG)
       zeta = (interpG - ig1_step) / (ig0 - ig1_step)
    end if
    ! returns interpG as zeta
    interpG = zeta

    if (inode == ionode .and. iprint_DM + min_layer >= 2) then
       do spin = 1, nspin
          write(io_lun, '(4x,a,i1," is: ",f16.6)') &
               trim(prefix)//" energy_1 for spin = ", spin, energy_1_spin(spin)
       end do
       write(io_lun, '(4x,a,f16.6)') trim(prefix)//" energy_1 overall is:      ", energy_1
    end if

    ! deallocate matrix
    do spin = nspin, 1, -1
       call free_temp_matrix(matM3old(spin))
    end do

    return
  end subroutine lineMinL
  !!***


  !!****f* DMMin/correct_electron_number
  !! PURPOSE
  !!   Selector subroutine for correct_electron_number
  !! INPUTS
  !! OUTPUT
  !! RETURN VALUE
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2011/08/11
  !! MODIFICATION HISTORY
  !!   2012/03/12 L.Tong
  !!   - rewrite for changed spin implementation
  !!   2015/06/08 lat
  !!    - Added experimental backtrace
  !! SOURCE
  !!
  subroutine correct_electron_number

    use mult_module,   only: matL, matphi
    use global_module, only: nspin, flag_fix_spin_population

    implicit none

    ! Local variables
    type(cq_timer) :: backtrace_timer
    integer        :: backtrace_level

    if (nspin == 1 .or. flag_fix_spin_population) then
       call correct_electron_number_fixspin
    else
       call correct_electron_number_varspin
    end if

    return
  end subroutine correct_electron_number
  !!*****


  !!****f* DMMin/correct_electron_number_fixspin
  !! PURPOSE
  !!   correct_electron_numbers for each spin component, treating
  !!   the spin populations fixed
  !! INPUTS
  !! OUTPUT
  !! RETURN VALUE
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2011/08/11
  !! MODIFICATION HISTORY
  !!   2012/03/11 L.Tong
  !!   - Major rewrite, the subroutine corrects electron number for
  !!     both spin channels.
  !!   2018/11/13 17:30 nakata
  !!    Changed matT and matTtran to be spin_SF dependent
  !! SOURCE
  !!
  subroutine correct_electron_number_fixspin

    use datatypes
    use logicals
    use numbers
    use matrix_data,    only: TLrange, Lrange
    use mult_module,    only: matT, matTtran, matL, matphi,           &
                              allocate_temp_matrix, free_temp_matrix, &
                              matrix_product_trace, matrix_product,   &
                              matrix_sum, mult, T_L_TL, TL_T_L,       &
                              LNV_matrix_multiply, matrix_transpose
    use primary_module, only: bundle
    use GenComms,       only: gsum, cq_abort
    use global_module,  only: ne_spin_in_cell, nspin, spin_factor, &
                              nspin_SF, flag_SpinDependentSF, min_layer
    use io_module,      only: return_prefix

    implicit none

    ! Local variables
    real(double)  :: step, step1, g0, g1, recA, B, C, D, truestep, dne
    integer       :: iter, spin, spin_SF
    logical       :: done
    real(double), dimension(nspin) :: electrons_0, electrons_2, &
                                      electrons, energy
    integer,      dimension(nspin) :: matTL, matphi2, matSphi, matSphi2
    character(len=20) :: subname = "correct_Ne: "
    character(len=120) :: prefix

    prefix = return_prefix(subname, min_layer)
    spin_SF = 1

    ! set electrons_0 to the correct electron number
    electrons_0(1:nspin) = ne_spin_in_cell(1:nspin)

    call matrix_transpose(matT(1), matTtran(1))
    if (nspin_SF==2) call matrix_transpose(matT(2), matTtran(2))
    do spin = 1, nspin
       done = .false.
       iter = 0
       if (flag_SpinDependentSF) spin_SF = spin
       ! allocate temporary matrices
       matTL(spin) = allocate_temp_matrix(TLrange,0)
       matphi2(spin) = allocate_temp_matrix(Lrange,0)
       matSphi(spin) = allocate_temp_matrix(Lrange,0)
       matSphi2(spin) = allocate_temp_matrix(Lrange,0)
       do while (.not. done .and. (iter < 20))
          ! get electron number and gradient
          call LNV_matrix_multiply(electrons, energy, dontK, dontM1,     &
                                   dontM2, dontM3, dontM4, dophi, dontE, &
                                   mat_phi=matphi, spin=spin)
          if (inode == ionode .and. iprint_DM + min_layer >= 2) &
               write(io_lun, '(4x,a,i1,") before correction: ",f16.6)') &
               trim(prefix)//" electron number (spin=", spin, electrons(spin)
          call matrix_product(matT(spin_SF), matphi(spin), matTL(spin), mult(T_L_TL))
          call matrix_product(matTL(spin), matTtran(spin_SF), matSphi(spin), mult(TL_T_L))
          g0 = matrix_product_trace(matSphi(spin), matphi(spin))
          ! initial guess is linear correction...
          step1 = (electrons_0(spin) - electrons(spin)) / g0
          step = 0.1_double
          if (inode == ionode .and. iprint_DM + min_layer >= 2) &
               write(io_lun, '(4x,a,i1,") are ",f16.6,e16.6)') &
                     trim(prefix)//" g0, step1 (spin=", spin, g0, step1
          if (abs(electrons_0(spin) - electrons(spin)) < 1.0e-9_double) then
             call matrix_sum(one, matL(spin), step1, matSphi(spin))
             ! check that electron number is correct
             if (iprint_DM + min_layer >= 2) then
                call LNV_matrix_multiply(electrons_2, energy, dontK,     &
                                         dontM1, dontM2, dontM3, dontM4, &
                                         dophi, dontE, mat_phi=matphi,   &
                                         spin=spin)
                if (inode == ionode) &
                     write(io_lun, '(4x,a,i1,") after correction: ",f15.6)') &
                                     trim(prefix)//" electron number (spin=", spin, electrons_2(spin)
             end if
             done = .true.
          else
             step = step1
             ! Take a step along the electron gradient
             call matrix_sum(one, matL(spin), step1, matSphi(spin))
             call LNV_matrix_multiply(electrons_2, energy, dontK, &
                                      dontM1, dontM2, dontM3, dontM4, &
                                      dophi, dontE, mat_phi=matphi2, &
                                      spin=spin)
             call matrix_product(matT(spin_SF), matphi2(spin), matTL(spin), &
                                 mult(T_L_TL))
             call matrix_product(matTL(spin), matTtran(spin_SF), &
                                 matSphi2(spin), mult(TL_T_L))
             g1 = matrix_product_trace(matphi(spin), matSphi2(spin))
             if (inode == ionode .and. iprint_DM + min_layer >= 2) &
                  write(io_lun, '(4x,a,i1,") are ",2f16.6)') &
                        trim(prefix)//" g1, elec2 (spin=", spin, g1, electrons_2(spin)
             ! get coefficients of polynomial
             D = electrons(spin) - electrons_0(spin)
             C = g0
             B = three * (electrons_2(spin) - electrons(spin)) / (step * step) - &
                 (g1 + two * g0) / step
             recA = (step * step * step) / &
                    (two * (electrons(spin) - electrons_2(spin)) + &
                     (g0 + g1) * step)
             B = B * recA
             C = C * recA
             D = D * recA
             truestep = SolveCubic(B, C, D, step, inode, ionode)
             if (inode .eq. ionode .and. iprint_DM + min_layer >= 2) &
                  write(io_lun, '(4x,a,i1,") are ",2e16.6)') &
                  trim(prefix)//" step, truestep (spin=", spin, step, truestep
             call matrix_sum(one, matL(spin), truestep - step, matSphi(spin))
             if (truestep == step) then
                if (inode == ionode .and. iprint_DM + min_layer >= 2) &
                     write(io_lun, '(4x,a)') trim(prefix)//" still in linear loop"
                ! check that electron number is correct
                call LNV_matrix_multiply(electrons_2, energy, dontK, &
                                         dontM1, dontM2, dontM3,     &
                                         dontM4, dophi, dontE,       &
                                         mat_phi=matphi, spin=spin)
                if (inode == ionode .and. iprint_DM + min_layer >= 2) &
                     write(io_lun,fmt='(4x,a,i1,") after correction",f16.6)') &
                     trim(prefix)//" electron number (spin=", spin, electrons_2(spin)
                dne = abs(electrons_2(spin) - electrons_0(spin))
                if (dne < 1e-4_double) done = .true.
                iter = iter + 1
             else
                done = .true.
                call LNV_matrix_multiply(electrons_2, energy, dontK, &
                                         dontM1, dontM2, dontM3,     &
                                         dontM4, dophi, dontE,       &
                                         mat_phi=matphi, spin=spin)
                if (inode == ionode .and. iprint_DM + min_layer >= 2) &
                     write(io_lun, '(4x,a,i1,") after correction:  ",f16.6)') &
                     trim(prefix)//" electron number (spin=", spin, electrons_2(spin)
             end if ! (truestep == step)
          end if ! (abs(electrons_0(spin) - electrons(spin)) < 1.0e-9_double)
       end do ! while ( .not. done .and. (iter < 20))
       call free_temp_matrix(matSphi2(spin))
       call free_temp_matrix(matSphi(spin))
       call free_temp_matrix(matphi2(spin))
       call free_temp_matrix(matTL(spin))
    end do ! spin

    return
  end subroutine correct_electron_number_fixspin
  !!*****


  !!****f* DMMin/correct_electron_number_varspin *
  !!
  !!  NAME
  !!   correct_electron_number_varspin
  !!  USAGE
  !!
  !!  PURPOSE
  !!   For spin polarised calculations:
  !!   Corrects the L matrices so that the electron number is correct
  !!   treating the overall L matrix as direct sum of the spin channels.
  !!
  !!   That means both spin channels share the same search steps, and
  !!   the steps and directions are calculated from the overall L
  !!   (direct sum of both spin components)
  !!
  !!   Note that the original number_of_bands in correct_electron_number
  !!   subroutine is an obsolete historical artifact. Use ne_in_cell
  !!   from the global module instead for the correct number of total
  !!   electrons.
  !!  INPUTS
  !!  USES
  !!  AUTHOR
  !!   L.Tong
  !!  CREATION DATE
  !!   2011/07/27
  !!  MODIFICATION HISTORY
  !!   Thursday, 2011/08/11 L.Tong
  !!     Fixed a bug, forgot to use matphi_dn from mult_module Switched
  !!     to LNV_matrix_multiply, as LNV_matrix_multiply_spin is now
  !!     obsolete
  !!   2018/11/13 17:30 nakata
  !!    Changed matT and matTtran to be spin_SF dependent
  !!  SOURCE
  !!
  subroutine correct_electron_number_varspin

    use datatypes
    use logicals
    use numbers
    use matrix_data,    only: TLrange, Lrange
    use mult_module,    only: matT, matTtran, matphi, matL,           &
                              allocate_temp_matrix, free_temp_matrix, &
                              matrix_product_trace, matrix_product,   &
                              matrix_sum, mult, T_L_TL, TL_T_L,       &
                              LNV_matrix_multiply, matrix_transpose
    use primary_module, only: bundle
    use GenComms,       only: gsum
    use global_module,  only: ne_in_cell, nspin, spin_factor, &
                              nspin_SF, flag_SpinDependentSF, min_layer
    use io_module,      only: return_prefix

    implicit none

    ! Local variables
    real(double), dimension(nspin) :: electrons_spin, energy_spin
    integer,      dimension(nspin) :: matTL, matphi2, matSphi, matSphi2
    real(double) :: electrons_0, electrons, energy, step, step1, &
                    electrons2, g0, g1, recA, B, C, D, truestep, dne
    integer :: iter, spin, spin_SF
    logical :: done
    character(len=20) :: subname = "correct_Ne: "
    character(len=120) :: prefix

    prefix = return_prefix(subname, min_layer)
    done = .false.
    iter = 0
    spin_SF = 1

    ! set the correct electron number
    electrons_0 = ne_in_cell

    do spin = 1, nspin
       matTL(spin) = allocate_temp_matrix(TLrange,0)
       matphi2(spin) = allocate_temp_matrix(Lrange,0)
       matSphi(spin) = allocate_temp_matrix(Lrange,0)
       matSphi2(spin) = allocate_temp_matrix(Lrange,0)
    end do

    ! get electron number and gradient
    call matrix_transpose(matT(1), matTtran(1))
    if (nspin_SF==2) call matrix_transpose(matT(2), matTtran(2))

    do while (.not. done .and. (iter < 20)) ! Was 20 !

       call LNV_matrix_multiply(electrons_spin, energy_spin, dontK, &
                                dontM1, dontM2, dontM3, dontM4, dophi,&
                                dontE, mat_phi=matphi)
       electrons = spin_factor * sum(electrons_spin)

       if (inode == ionode .and. iprint_DM + min_layer >= 2) &
            write(io_lun, '(4x,a,3f16.6)') &
            trim(prefix)//" electron numbers before correction (up, down, total): ", &
            electrons_spin(1), electrons_spin(nspin), electrons
       g0 = zero
       do spin = 1, nspin
          if (flag_SpinDependentSF) spin_SF = spin
          call matrix_product(matT(spin_SF), matphi(spin), matTL(spin), mult(T_L_TL))
          call matrix_product(matTL(spin), matTtran(spin_SF), matSphi(spin), mult(TL_T_L))
          g0 = g0 + spin_factor * &
               matrix_product_trace(matSphi(spin), matphi(spin))
       end do

       ! initial guess is linear correction...
       step1 = (electrons_0 - electrons) / (g0 + RD_ERR)
       step = 0.1_double
       if (inode == ionode .and. iprint_DM + min_layer >= 2) &
            write(io_lun, '(4x,a,f16.6,e16.6)') trim(prefix)//" g0, step1 are", g0, step1

       ! if we are within 0.1% of the correct number, linear will do.
       if (abs(electrons_0 - electrons) < 1.0e-9_double) then
          do spin = 1, nspin
             call matrix_sum(one, matL(spin), step1, matSphi(spin))
          end do
          done = .true.
       else
          step = step1
          ! Take a step along the electron gradient
          do spin = 1, nspin
             call matrix_sum (one, matL(spin), step1, matSphi(spin))
          end do

          call LNV_matrix_multiply(electrons_spin, energy_spin, dontK,&
                                   dontM1, dontM2, dontM3, dontM4,    &
                                   dophi, dontE, mat_phi=matphi2)
          electrons2 = spin_factor * sum(electrons_spin)

          g1 = zero
          do spin = 1, nspin
             if (flag_SpinDependentSF) spin_SF = spin
             call matrix_product(matT(spin_SF), matphi2(spin), matTL(spin), &
                                 mult(T_L_TL))
             call matrix_product(matTL(spin), matTtran(spin_SF), &
                                 matSphi2(spin), mult(TL_T_L))
             g1 = g1 + spin_factor * &
                  matrix_product_trace(matphi(spin), matSphi2(spin))
          end do
          if (inode == ionode .and. iprint_DM + min_layer >= 2)               &
               write(io_lun, '(4x,a,4f16.6)') &
               trim(prefix)//" g1, elec (up, down, total): ", &
               g1, electrons_spin(1), electrons_spin(nspin),   &
               electrons2

          ! get coefficients of polynomial
          D = electrons - electrons_0
          C = g0
          !ORI B = (3.0_double*electrons2 - g1*step - &
          !ORI     2.0_double*g0*step - 3.0_double*electrons) / &
          !ORI     (step*step)
          B = three * (electrons2 - electrons) / (step * step) - &
               (g1 + two * g0) / step
          !ORI recA = (step*step*step) / &
          !ORI      (2.0_double*electrons - 2.0_double*electrons2 + &
          !ORI      g0*step + g1*step)
          recA = (step * step * step) / &
               (two * (electrons - electrons2) + (g0 + g1) * step)
          B = B * recA
          C = C * recA
          D = D * recA

          truestep = SolveCubic(B, C, D, step, inode, ionode)
          if (inode == ionode .and. iprint_DM + min_layer >= 2) &
               write(io_lun, '(4x,a,2e16.6)') trim(prefix)//" step, truestep ", step, truestep

          do spin = 1, nspin
             call matrix_sum(one, matL(spin), truestep - step, matSphi(spin))
          end do

          if (truestep == step) then

             if (inode == ionode .and. iprint_DM + min_layer >= 2) &
                  write(io_lun, '(4x,a)') trim(prefix)//" still in linear loop"

             ! check that electron number is correct
             call LNV_matrix_multiply(electrons_spin, energy_spin,   &
                                      dontK, dontM1, dontM2, dontM3, &
                                      dontM4, dophi, dontE, mat_phi=matphi)
             electrons2 = spin_factor * sum(electrons_spin)

             if (inode == ionode .and. iprint_DM + min_layer >= 2) &
                  write(io_lun,fmt='(4x,a,3f16.6)') &
                  trim(prefix)//" electron numbers after correction (up, down, total):  ", &
                  electrons_spin(1), electrons_spin(nspin), electrons2
             dne = abs(electrons2 - electrons_0)
             if ((dne / electrons_0) < 1.0e-6_double) done = .true.
             iter = iter + 1
          else
             done = .true.
             call LNV_matrix_multiply(electrons_spin, energy_spin, &
                                      dontK, dontM1, dontM2, dontM3, &
                                      dontM4, dophi, dontE, mat_phi=matphi)
             electrons2 = spin_factor * sum(electrons_spin)
             if (inode == ionode .and. iprint_DM + min_layer >= 2) &
                  write(io_lun,fmt='(4x,a,3f16.6)') &
                  trim(prefix)//" electron numbers after correction (up, down, total):  ", &
                  electrons_spin(1), electrons_spin(nspin), electrons2
          end if
       end if ! (abs(electrons_0 - electrons) < 1.0e-9_double)
    end do ! while (.not. done .and. iter < 20)

    ! free temp matrices
    do spin = nspin, 1, -1
       call free_temp_matrix (matSphi2(spin))
       call free_temp_matrix (matSphi(spin))
       call free_temp_matrix (matphi2(spin))
       call free_temp_matrix (matTL(spin))
    end do

    return

1   format(2x,'Electron numbers before correction (Nup, Ndn, Ntotal) :', &
           f25.15, f25.15, f25.15)
2   format(2x,'Electron number after correction (Nup, Ndn, Ntotal) :', &
           f25.15, f25.15, f25.15)

  end subroutine correct_electron_number_varspin
  !!****

! -----------------------------------------------------------
! Function SolveCubic
! -----------------------------------------------------------

!!****f* DMMin/SolveCubic *
!!
!!  NAME
!!   SolveCubic
!!  USAGE
!!
!!  PURPOSE
!!   Solves a cubic: x^3+Ax^2+Bx+C
!!   This is a general solution first worked out in the
!!   sixteenth century (published by Gerolamo Cardano)
!!
!!   See, for example, Abramowitz & Stegun p. 17
!!  INPUTS
!!
!!
!!  USES
!!
!!  AUTHOR
!!   E.H.Hernandez
!!  CREATION DATE
!!   10/04/95 ?
!!  MODIFICATION HISTORY
!!   30/05/2001 dave
!!    Added ROBODoc header
!!   20/06/2001 dave
!!    Included in DMMinModule
!!   2019/10/23 16:57 dave
!!    Updated and polished
!!  SOURCE
!!
  function SolveCubic(a,b,c,guess,inode,ionode)

    use datatypes
    use numbers
    use global_module, only: min_layer

    implicit none

    integer inode,ionode
    real(double) :: SolveCubic
    real(double) :: a,b,c,guess

    real(double) :: Q, R, S, T, theta, z1, z2, z3

    Q = (a*a - three*b) / nine
    R = (two*a*a*a - nine*a*b + 27.0_double*c)/54.0_double

    if ((R*R)<(Q*Q*Q)) then ! Guarantees Q is positive
       ! three roots...
       theta = acos(R/(sqrt(Q*Q*Q)))
       T = -two * sqrt(Q)
       ! Note that solutions are given relative to guess (step taken to find parameters)
       z1 = T * cos(theta/three) - (a/three) 
       z2 = T * cos((theta+twopi)/three) - (a/three) 
       z3 = T * cos((theta-twopi)/three) - (a/three) 

       ! Take solution with smallest magnitude
       if (abs(z1)<=abs(z2).and.abs(z1)<=abs(Z3)) then
          SolveCubic = z1 
       else if (abs(z2)<=abs(z3)) then
          SolveCubic = z2 
       else
          SolveCubic = z3 
       end if
    else

       ! the solution for the single real root
       Q = -Q
       R = -R

       S = R + sqrt(Q*Q*Q + R*R)
       T = R - sqrt(Q*Q*Q + R*R)

       if ( S >= zero ) then
          S = S**third
       else
          S = -abs(S)**third
       end if

       if ( T >= zero ) then
          T = T**third
       else
          T = -abs(T)**third
       end if

       SolveCubic = S + T - a/three
       ! Experience suggests that, while this is a correct solution to the cubic,
       ! a large step is completely wrong for electron number correction - test
       ! for this, and if necessary return only the initial guess
       ! 
       ! The prefactor of 10 is arbitrary but seems reasonable
       if(abs(SolveCubic)>10.0_double*abs(guess)) then
          if(inode==ionode .and. iprint_DM + min_layer > 2) &
               write (io_lun, '(8x,"Step too large for Ne correction: linear guess was: ",2e16.6)') &
                     SolveCubic, guess
          SolveCubic = guess
       end if

    end if

    return
  end function SolveCubic
!!***
end module DMMin
