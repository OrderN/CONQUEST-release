! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
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
!!   Contains initialisation, early stage (steepest descents) and
!!   late stage (GR-Pulay - D.R.Bowler and M.J.Gillan, Chem. Phys.
!!   Lett. 325, 473 (2000)) techniques as well as a call to a diagonalisation routine for exact solution
!!  USES
!!   datatypes, DiagModule, GenBlas, GenComms, global_module, logicals, matrix_data, maxima_module, McWeeny,
!!   multiply_module, mult_module, numbers, PosTan, primary_module, Pulay
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
!!    Minor changes to shift H back after DM minimisation, and to switch between O(N) and diagonalisation
!!   14:23, 26/02/2003 drb 
!!    Added electron number check to linear part of correct_electron_number
!!   16:43, 10/05/2005 dave 
!!    Various small changes throughout code: mainly bug fixes
!!   2008/02/01 17:46 dave
!!    Changes to write output to file not stdout
!!   2008/05/22 ast
!!    Added timers
!!   2012/03/01 L.Tong
!!    Added interface for lineMinL
!!  SOURCE
!!
module DMMin

  use datatypes
  use global_module, ONLY: io_lun
  use timer_stdclocks_module, ONLY: start_timer,stop_timer,tmr_std_densitymat

  integer :: maxpulayDMM
  integer,save :: n_dumpL = 5
  real(double) :: LinTol_DMM
  !integer, parameter :: mx_pulay = 5
  !real(double), parameter :: LinTol = 0.1_double

  !!****f* DMMin/lineMinL
  !! PURPOSE
  !!   Interface for lineMinL_nospin and lineMinL_spin
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/03/01
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  interface lineMinL
     module procedure lineMinL_nospin
     module procedure lineMinL_spin
  end interface lineMinL
  !!*****
  
  ! -------------------------------------------------------
  ! RCS ident string for object file id
  ! -------------------------------------------------------
  character(len=80), save, private :: RCSid = "$Id$"
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
!!  SOURCE
!!
  subroutine FindMinDM (n_L_iterations, vary_mu, tolerance, mu,&
       & inode, ionode, resetL, record)

    use datatypes
    use numbers
    use global_module, ONLY: iprint_DM, IPRINT_TIME_THRES1, &
         flag_spin_polarisation, flag_fix_spin_population, ne_in_cell,&
         ne_up_in_cell, ne_dn_in_cell
    use mult_module, ONLY: matrix_transpose, matT, matTtran, matL, &
         matL_dn, matrix_sum
    use McWeeny, ONLY: InitMcW, McWMin
    use PosTan, ONLY: max_iters, cscale, PulayE, PulayR, PulayC, &
         PulayBeta, pos_tan, fit_coeff
    use DiagModule, ONLY: FindEvals, diagon
    use io_module, ONLY: dump_matrix
    use energy, ONLY: entropy
    use timer_module, ONLY: cq_timer,start_timer,stop_print_timer,&
         WITH_LEVEL

    implicit none

    ! Passed variables
    integer :: n_L_iterations, inode, ionode
    real(double) :: mu, tolerance
    logical :: resetL, vary_mu, record

    ! Local variables
    logical :: done, inflex, problem, early
    integer :: ndone, nkeep, ndelta, i
    real(double) :: LastE
    real(double), parameter :: eprec = 1.0e-12_double
    real(double), save :: delta_e
    type(cq_timer) :: tmr_l_iter
!TM 2010.Nov.06
    integer :: niter = 0
    integer, parameter :: niter_max=10

    call start_timer(tmr_std_densitymat)
    entropy = zero
    if(diagon) then ! Use exact diagonalisation to get K
       if (flag_fix_spin_population) then
          call FindEvals(ne_in_cell, ne_up_in_cell, ne_dn_in_cell)
       else
          call FindEvals(ne_in_cell)
       end if
       call stop_timer(tmr_std_densitymat)
       return
    end if
    if(inode==ionode) &
         write(io_lun,*) 'Welcome to FindDMMin, tol: ', &
         tolerance, n_L_iterations
    problem = .false.
    inflex = .false.
    done = .false.
    early = .false.
    ! Start with the transpose of S^-1
    call matrix_transpose(matT,matTtran)
    ! Now minimise the energy
    ndone = 0
    niter = 0
    do while(.NOT.done)
       call start_timer(tmr_l_iter,WITH_LEVEL)
       if(resetL.OR.inflex) then ! Reset L to McW
          call InitMcW (inode, ionode)
          call McWMin (n_L_iterations, delta_e) 
          early = .true.
          resetL = .false.    !2010.Nov.06 TM
       endif
       if(early.OR.problem) then
          inflex = .false.
          call earlyDM (ndone, n_L_iterations, delta_e, done,&
               & vary_mu, tolerance, mu, inode, ionode, inflex,&
               & record)
       endif
       if((.NOT.done).AND.(.NOT.inflex)) then ! Continue on to late stage
          call lateDM (ndone, n_L_iterations, done, delta_e, vary_mu,&
               & mu, inode, ionode, tolerance, record)
       endif
       niter = niter + 1
       if(problem) resetL = .true.
       !ORI problem = .true. ! If we're not done, then this kicks us back to early
       if(niter > niter_max) then
         problem = .true. ! If we're not done, then this kicks us back to early
         niter = 0
       endif
       call stop_print_timer(tmr_l_iter,"FindMinDM iteration",IPRINT_TIME_THRES1)
    enddo ! end of do while (.NOT.done)
    ! *** Add frequency of output here *** !

    call dump_matrix("L",matL,inode)
    if (flag_spin_polarisation) call dump_matrix ("L_dn", matL_dn, inode)
 
    if(record) then
       if(inode==ionode.AND.iprint_DM>1) then
          write(io_lun,*) '  List of residuals and energies'
          do i=1,ndone
             write(io_lun,7) i, PulayR(i), PulayE(i)
          enddo
       endif
       call fit_coeff(PulayC, PulayBeta,PulayE, PulayR,ndone)
       if(inode==ionode) write(io_lun,6) PulayC,PulayBeta
    endif
    call stop_timer(tmr_std_densitymat)
6   format(2x,'dE to dR parameters - C: ',f15.8,' beta: ',f15.8)
7   format(2x,i4,2f15.8)
8   format(2x,'Welcome to minimise_energy. Tolerance is ',f15.8)
    return
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
!!    Added call to move H back (i.e. undo the H-mu.S shift) after minimisation
!!   2004/10/28 drb
!!    Removed shift of H
!!   10:09, 13/02/2006 drb 
!!    Removed all explicit references to data_ variables and rewrote in terms of new 
!!    matrix routines
!!   2009/07/08 16:41 dave
!!    Introduced atom-based tolerance (just divide by ni_in_cell)
!!   2011/08/14 L.Tong
!!    Added spin polarisation
!!    Added matphi_dn
!!   2011/08/25 L.Tong
!!    Removed the redundant parameter number_of_bands
!!  SOURCE
!!
  subroutine earlyDM (ndone, n_L_iterations, delta_e, done, vary_mu,&
       & tolerance, mu, inode, ionode, inflex, record)

    use datatypes
    use logicals
    use numbers
    use matrix_data, ONLY: Lrange, TLrange
    use mult_module, ONLY: matT, matphi, matphi_dn, mult, T_L_TL,&
         & TL_T_L, matrix_product_trace, LNV_matrix_multiply,&
         & allocate_temp_matrix, free_temp_matrix, matrix_sum,&
         & matrix_product
    use primary_module, ONLY: bundle
    use PosTan, ONLY: PulayR, PulayE
    use GenComms, ONLY: cq_abort, gsum
    use global_module, ONLY: iprint_DM,IPRINT_TIME_THRES1, ni_in_cell,&
         flag_global_tolerance, flag_spin_polarisation, &
         flag_fix_spin_population
    use timer_module, ONLY: cq_timer, start_timer, stop_print_timer,&
         & WITH_LEVEL

    implicit none

    ! Passed variables
    real(double) :: mu, tolerance, delta_e
    integer :: inode, ionode, n_L_iterations, ndone
    logical :: inflex, vary_mu, done, record

    ! Local variables
    integer :: n_iter, length
    integer :: matM3, matSM3, matSphi, mat_temp, mat_search
    integer :: matM3_dn, matSM3_dn, mat_search_dn, matSphi_dn
    real(double) :: energy0, energy0_up, energy0_dn 
    real(double) :: energy1, energy1_up, energy1_dn
    real(double) :: electrons, electrons_up, electrons_dn
    real(double) :: e_dot_n, e_dot_n_up, e_dot_n_dn
    real(double) :: n_dot_n, n_dot_n_up, n_dot_n_dn
    real(double) :: g0, g1
    real(double) ::  interpG,  zeta, direct_sum_factor
    type(cq_timer) :: tmr_l_tmp1, tmr_l_iter

    call start_timer(tmr_l_tmp1,WITH_LEVEL)
    matM3 = allocate_temp_matrix(Lrange,0)
    matSM3 = allocate_temp_matrix(Lrange,0)
    matSphi = allocate_temp_matrix(Lrange,0)
    mat_temp = allocate_temp_matrix(TLrange,0)
    mat_search = allocate_temp_matrix(Lrange,0)
    if (flag_spin_polarisation) then
       matM3_dn = allocate_temp_matrix (Lrange,0)
       matSM3_dn = allocate_temp_matrix (Lrange,0)
       matSphi_dn = allocate_temp_matrix (Lrange,0)
       mat_search_dn = allocate_temp_matrix (Lrange,0)
    end if
    if(ndone>n_L_iterations) call cq_abort('earlyDM: too many L&
         & iterations', ndone, n_L_iterations)
    ! Set up the gradient and electron number
    if (vary_mu) then
       call correct_electron_number (iprint_DM, inode, ionode)
    endif    
    if (flag_spin_polarisation) then 
       call LNV_matrix_multiply (electrons_up, energy0_up, dontK, &
            dontM1, dontM2, doM3, dontM4, dophi, doE, 0, matM3, 0, &
            matphi, spin=1)
       call LNV_matrix_multiply (electrons_dn, energy0_dn, dontK, &
            dontM1, dontM2, doM3, dontM4, dophi, doE, 0, matM3_dn, 0, &
            matphi_dn, spin=2)
       electrons = electrons_up + electrons_dn
       ! energy0 = E(L_n_iter)
       energy0 = energy0_up + energy0_dn
       ! Pre- and post-multiply M3 by S^-1 so that it is contravariant
       call matrix_product (matT, matM3, mat_temp, mult(T_L_TL))
       call matrix_product (mat_temp, matT, matSM3, mult(TL_T_L))
       call matrix_product (matT, matM3_dn, mat_temp, mult(T_L_TL))
       call matrix_product (mat_temp, matT, matSM3_dn, mult(TL_T_L))
       ! Project gradient perpendicular to electron gradient
       call matrix_sum (zero, mat_search, -one, matSM3)
       call matrix_sum (zero, mat_search_dn, -one, matSM3_dn)
       if (vary_mu) then
          ! Pre- and post-multiply phi by S^-1 so that it is contravariant
          call matrix_product (matT, matphi, mat_temp, mult(T_L_TL))
          call matrix_product (mat_temp, matT, matSphi, mult(TL_T_L))
          call matrix_product (matT, matphi_dn, mat_temp, mult(T_L_TL))
          call matrix_product (mat_temp, matT, matSphi_dn, mult(TL_T_L))
          ! Only one is pre- and post-multiplied because (A,B) = Tr(ASBS)
          e_dot_n_up = matrix_product_trace (matSM3, matphi)
          e_dot_n_dn = matrix_product_trace (matSM3_dn, matphi_dn)
          e_dot_n = e_dot_n_up + e_dot_n_dn
          n_dot_n_up = matrix_product_trace (matSphi, matphi)
          n_dot_n_dn = matrix_product_trace (matSphi_dn, matphi_dn)
          n_dot_n = n_dot_n_up + n_dot_n_dn
          if (inode .eq. ionode .and. iprint_DM >= 2) then
             write (io_lun, *) 'e_up.n_up and e_dn.n_dn', e_dot_n_up, e_dot_n_dn 
             write (io_lun, *) 'n_up.n_up and n_dn.n_dn', n_dot_n_up, n_dot_n_dn 
          end if
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
          if (flag_fix_spin_population) then
             call matrix_sum (one, mat_search, (e_dot_n_up) / &
                  (n_dot_n_up), matSphi)
             call matrix_sum (one, mat_search_dn, (e_dot_n_dn) /&
                  (n_dot_n_dn), matSphi_dn)
          else
             direct_sum_factor = e_dot_n / n_dot_n
             call matrix_sum (one, mat_search, direct_sum_factor, matSphi)
             call matrix_sum (one, mat_search_dn, direct_sum_factor, matSphi_dn)
          end if
       endif ! if (vary_mu) then
    else ! if (flag_spin_ploarisation) then
       call LNV_matrix_multiply (electrons, energy0, dontK, dontM1, &
            dontM2, doM3, dontM4, dophi, doE, 0, matM3, 0, matphi)
       ! Pre- and post-multiply M3 by S^-1 so that it is contravariant
       call matrix_product(matT, matM3, mat_temp, mult(T_L_TL))
       call matrix_product(mat_temp, matT, matSM3, mult(TL_T_L))
       ! Project gradient perpendicular to electron gradient
       call matrix_sum(zero, mat_search, -one, matSM3)
       if(vary_mu) then
          ! Pre- and post-multiply phi by S^-1 so that it is contravariant
          call matrix_product(matT, matphi, mat_temp, mult(T_L_TL))
          call matrix_product(mat_temp, matT, matSphi, mult(TL_T_L))
          ! Only one is pre- and post-multiplied because (A,B) = Tr(ASBS)
          e_dot_n = matrix_product_trace(matSM3, matphi)
          n_dot_n = matrix_product_trace(matSphi, matphi)
          !mu = e_dot_n/n_dot_n
          if(inode.eq.ionode.and.iprint_DM>=2) &
               write(io_lun,*) 'e.n and n.n: ', e_dot_n, n_dot_n
          ! This is right: search dirn is -SM3
          call matrix_sum(one,mat_search, (e_dot_n)/(n_dot_n), matSphi)
          ! Here, we can't alter M3 directly: lineMinL expects REAL gradient
       endif ! if (vary_mu) then
    end if ! if (flag_spin_ploarisation) then
    call stop_print_timer(tmr_l_tmp1, "earlyDM - Preliminaries", &
         IPRINT_TIME_THRES1)
    !-----------!
    ! MAIN LOOP !
    !-----------!
    ! delta_e is used for spin polarised calculations too
    delta_e = zero
    if (flag_spin_polarisation) then
       do n_iter = 1, n_L_iterations
          call start_timer (tmr_l_iter,WITH_LEVEL)
          ! g0 = tr(matM3(L_n_iter) * matSM3(L_n_iter))
          if (flag_global_tolerance) then
             g0 = matrix_product_trace (matM3, matSM3)
             g0 = g0 + matrix_product_trace (matM3_dn, matSM3_dn)
          else
             g0 = matrix_product_trace (matM3, matSM3) / &
                  real (ni_in_cell, double)
             g0 = g0 + matrix_product_trace (matM3_dn, matSM3_dn) / &
                  real (ni_in_cell, double)
          end if
          if (inode == ionode .and. iprint_DM >= 1) &
               write (io_lun, 1) n_iter, energy0, g0
          ! minimise total E along search direction, updates energy1 =
          ! E(L_n_iter+1), delta_e = energy1 - energy0,
          ! matM3(L_n_iter+1), matM3_dn(L_n_iter+1), inflex and
          ! interpG
          call lineMinL(iprint_DM, matM3, matM3_dn, mat_search, &
               mat_search_dn, mat_temp, matSM3, matSM3_dn, energy0, &
               energy1, delta_e, inode, ionode, inflex, interpG)
          ! panic if found inflexion point
          if (inflex) then
             if (inode == ionode) write (io_lun,*) 'Panic ! Inflexion point found !'
             ndone = n_iter
             call free_temp_matrix (mat_search_dn)
             call free_temp_matrix (matSphi_dn)
             call free_temp_matrix (matSM3_dn)
             call free_temp_matrix (matM3_dn)
             call free_temp_matrix (mat_search)
             call free_temp_matrix (mat_temp)
             call free_temp_matrix (matSphi)
             call free_temp_matrix (matSM3)
             call free_temp_matrix (matM3)
             call stop_print_timer (tmr_l_iter,"an earlyDM iteration"&
                  &,IPRINT_TIME_THRES1)
             return !This is a panic sign !
          endif
          ! Gradient after line min - assumes LVN_matrix_mult at end
          !  of lineMinL
          ! Pre- and post-multiply M3 by S^-1 so that it is
          !  contravariant, matSM3 = matSM3(L_n_iter+1)
          call matrix_product (matT, matM3, mat_temp, mult(T_L_TL))
          call matrix_product (mat_temp, matT, matSM3, mult(TL_T_L))
          call matrix_product (matT, matM3_dn, mat_temp, mult(T_L_TL))
          call matrix_product (mat_temp, matT, matSM3_dn, mult(TL_T_L))
          ! g1 = tr(matM3(L_n_iter+1) * matSM3(L_n_iter+1))
          if(flag_global_tolerance) then
             g1 = matrix_product_trace (matM3, matSM3)
             g1 = g1+ matrix_product_trace (matM3_dn, matSM3_dn)
          else
             g1 = matrix_product_trace (matM3, matSM3) / &
                  real (ni_in_cell,double)
             g1 = g1 + matrix_product_trace (matM3_dn, matSM3_dn) / &
                  real (ni_in_cell,double)
          end if
          ! Test for linearity or convergence
          ! zeta is returned by lineMinL
          zeta = interpG
          if (inode == ionode .and. iprint_DM >= 2)&
               & write (io_lun, 2) delta_e, zeta
          if (inode == ionode .and. iprint_DM >= 2)&
               & write (io_lun, *) 'zeta...: ', g0, g1, interpG
          ! Correct L_n_iter+1 again to get electron numbers correct
          if (vary_mu) then
             call correct_electron_number (iprint_DM, inode, ionode)
          endif
          ! recalculate the quantities after L_n_iter+1 is updated
          ! 2011/08/24 L.Tong: 
          !   The orignal implementation uses energy0 as input for
          !   LNV_matrix_multiply, but I think this is incorrect and
          !   should be energy1
          call LNV_matrix_multiply (electrons_up, energy1_up, dontK, dontM1,&
               & dontM2, doM3, dontM4, dophi, doE, 0, matM3, 0, matphi, spin=1)
          call LNV_matrix_multiply (electrons_dn, energy1_dn, dontK,&
               & dontM1, dontM2, doM3, dontM4, dophi, doE, 0,&
               & matM3_dn, 0, matphi_dn, spin=2)
          electrons = electrons_up + electrons_dn
          energy1 = energy1_up + energy1_dn
          ! Pre- and post-multiply M3 by S^-1 so that it is
          ! contravariant
          call matrix_product (matT, matM3, mat_temp, mult(T_L_TL))
          call matrix_product (mat_temp, matT, matSM3, mult(TL_T_L))
          call matrix_product (matT, matM3_dn, mat_temp, mult(T_L_TL))
          call matrix_product (mat_temp, matT, matSM3_dn, mult(TL_T_L))
          ! prepare for the next iterative step
          energy0 = energy1
          ! update the search direction
          call matrix_sum (zero, mat_search, -one, matSM3)
          call matrix_sum (zero, mat_search_dn, -one, matSM3_dn)
          ! Project search direction perpendicular to electron gradient
          if (vary_mu) then
             ! Wrap electron gradient for tensorial correctness
             call matrix_product (matT, matphi, mat_temp, mult(T_L_TL))
             call matrix_product (mat_temp, matT, matSphi, mult(TL_T_L))
             call matrix_product (matT, matphi_dn, mat_temp, mult(T_L_TL))
             call matrix_product (mat_temp, matT, matSphi_dn, mult(TL_T_L))
             e_dot_n_up = matrix_product_trace (matSM3, matphi)
             e_dot_n_dn = matrix_product_trace (matSM3_dn, matphi_dn)
             e_dot_n = e_dot_n_up + e_dot_n_dn
             n_dot_n_up = matrix_product_trace (matSphi, matphi)
             n_dot_n_dn = matrix_product_trace (matSphi_dn, matphi_dn)
             n_dot_n = n_dot_n_up + n_dot_n_dn
             if (flag_fix_spin_population) then
                call matrix_sum (one, mat_search, &
                     (e_dot_n_up) / (n_dot_n_up), matSphi)
                call matrix_sum (one, mat_search_dn, &
                     (e_dot_n_dn) /(n_dot_n_dn), matSphi_dn)
             else
                direct_sum_factor = e_dot_n / n_dot_n
                call matrix_sum (one, mat_search, direct_sum_factor, matSphi)
                call matrix_sum (one, mat_search_dn, direct_sum_factor, matSphi_dn)
             end if
             ! Here, we can't alter M3 directly: lineMinL expects REAL gradient
          end if
          if (flag_global_tolerance) then
             g0 = matrix_product_trace (matM3, matSM3)
             g0 = g0 + matrix_product_trace (matM3_dn, matSM3_dn)
          else
             g0 = matrix_product_trace (matM3, matSM3) / &
                  real (ni_in_cell, double)
             g0 = g0 + matrix_product_trace (matM3_dn, matSM3_dn) / &
                  real (ni_in_cell, double)
          end if
          ! store the Pulay histories
          PulayR(ndone + n_iter) = g0
          PulayE(ndone + n_iter) = energy0
          ! test if we are linear
          if (abs (zeta) < LinTol_DMM) then ! We're linear Added +1 to
             ! n_iter for cosmetic reasons (this is true since at this
             ! stage g0 and energy0 etc are already prepared for
             ! n_iter+1 step)
             if (inode == ionode .and. iprint_DM >= 1)&
                  & write (io_lun, 1) n_iter+1, energy0, g0
             if ((inode == ionode).and. (iprint_DM >= 2))&
                  & write (io_lun, *) 'Linearity satisfied - calling PulayL'
             ndone = n_iter
             call free_temp_matrix (mat_search_dn)
             call free_temp_matrix (matSphi_dn)
             call free_temp_matrix (matSM3_dn)
             call free_temp_matrix (matM3_dn)
             call free_temp_matrix (mat_search)
             call free_temp_matrix (mat_temp)
             call free_temp_matrix (matSphi)
             call free_temp_matrix (matSM3)
             call free_temp_matrix (matM3)
             call stop_print_timer (tmr_l_iter, "an earlyDM&
                  & iteration", IPRINT_TIME_THRES1)
             return
          endif
          call stop_print_timer (tmr_l_iter, "an earlyDM iteration"&
               &, IPRINT_TIME_THRES1)
       end do ! do n_iter = 1, n_L_iterations
    else ! if (flag_spin_polarisation) then
       do n_iter = 1, n_L_iterations
          call start_timer(tmr_l_iter,WITH_LEVEL)
          ! Gradient before line min
          if(flag_global_tolerance) then
             g0 = matrix_product_trace(matM3, matSM3)
          else
             g0 = matrix_product_trace(matM3, matSM3) / &
                  real(ni_in_cell, double)
          end if
          if (inode == ionode .AND. iprint_DM >= 1) &
               write (io_lun, 1) n_iter, energy0, g0
          call lineMinL(iprint_DM, matM3, mat_search, mat_temp, &
               matSM3, energy0, energy1, delta_e, inode, ionode, &
               inflex, interpG)
          ! 2011/08/17 LT: delta_e is already calculated in lineMinL,
          ! the following line redundant?
          ! delta_e = energy1 - energy0
          if (inflex) then
             if (inode == ionode) &
                  write (io_lun, *) 'Panic ! Inflexion point found !'
             ndone = n_iter
             call free_temp_matrix(mat_search)
             call free_temp_matrix(mat_temp)
             call free_temp_matrix(matSphi)
             call free_temp_matrix(matSM3)
             call free_temp_matrix(matM3)
             call stop_print_timer(tmr_l_iter, "an earlyDM iteration",&
                  IPRINT_TIME_THRES1)
             return !This is a panic sign !
          endif
          ! Gradient after line min - assumes main_matrix_mult at end of lineMinL
          ! Pre- and post-multiply M3 by S^-1 so that it is contravariant
          call matrix_product(matT, matM3, mat_temp, mult(T_L_TL))
          call matrix_product(mat_temp, matT, matSM3, mult(TL_T_L))       
          if (flag_global_tolerance) then
             g1 = matrix_product_trace(matM3, matSM3)
          else
             g1 = matrix_product_trace(matM3, matSM3) / &
                  real(ni_in_cell, double)
          end if
          ! Test for linearity or convergence
          !zeta = (interpG - g1)/(g0)
          ! DRB 2004/09/28 Now zeta returned by lineMinL
          zeta = interpG
          if (inode == ionode .AND. iprint_DM >= 2) &
               write (io_lun, 2) delta_e, zeta
          if (inode == ionode .AND. iprint_DM >= 2) &
               write (io_lun, *) 'zeta...: ', g0, g1, interpG
          if (vary_mu) then
             call correct_electron_number(iprint_DM, inode, ionode)
          endif
          ! 2011/08/24 L.Tong: 
          !   The orignal implementation uses energy0 as input for
          !   LNV_matrix_multiply, but I think this is incorrect and
          !   should be energy1
          call LNV_matrix_multiply(electrons, energy1, dontK, dontM1, &
               dontM2, doM3, dontM4, dophi, doE, 0, matM3, 0, matphi)
          ! Pre- and post-multiply M3 by S^-1 so that it is
          !  contravariant
          call matrix_product(matT, matM3, mat_temp, mult(T_L_TL))
          call matrix_product(mat_temp, matT, matSM3, mult(TL_T_L))
          energy0 = energy1
          call matrix_sum(zero, mat_search, -one, matSM3)
          ! Project search direction perpendicular to electron gradient
          if (vary_mu) then
             ! Wrap electron gradient for tensorial correctness
             call matrix_product(matT, matphi, mat_temp, mult(T_L_TL))
             call matrix_product(mat_temp, matT, matSphi, mult(TL_T_L))
             e_dot_n = matrix_product_trace(matSM3, matphi)
             n_dot_n = matrix_product_trace(matSphi, matphi)
             call matrix_sum(one, mat_search, (e_dot_n) / (n_dot_n), matSphi)
             ! Here, we can't alter M3 directly: lineMinL expects REAL gradient
          end if
          if (flag_global_tolerance) then
             g0 = matrix_product_trace(matM3, matSM3)
          else
             g0 = matrix_product_trace(matM3, matSM3) / &
                  real(ni_in_cell, double)
          end if
          PulayR(ndone + n_iter) = g0
          PulayE(ndone + n_iter) = energy0
          if (abs(zeta) < LinTol_DMM) then ! We're linear
             ! Commented out test for convergence to force checking of 
             ! convergence by Pulay - DRB, 22/03/01
             ! Added +1 to n_iter for cosmetic reasons
             if (inode == ionode .AND. iprint_DM >= 1) &
                  write (io_lun, 1) n_iter+1, energy0, g0
             if ((inode == ionode) .AND. (iprint_DM >= 2)) &
                  write (io_lun, *) 'Linearity satisfied - calling PulayL'
             !if(abs(delta_E)<tolerance) then
             !  if((inode==ionode).AND.(iprint_DM>=2)) &
             !    write(io_lun,*) 'Tolerance also achieved'
             !  done = .true.
             !endif
             ndone = n_iter
             call free_temp_matrix(mat_search)
             call free_temp_matrix(mat_temp)
             call free_temp_matrix(matSphi)
             call free_temp_matrix(matSM3)
             call free_temp_matrix(matM3)
             call stop_print_timer(tmr_l_iter, "an earlyDM iteration",&
                  IPRINT_TIME_THRES1)
             return
             !else if(abs(delta_E)<tolerance) then
             !   if(inode==ionode.AND.iprint_DM>=1) write(io_lun,1) n_iter, energy0, g0
             !   if((inode==ionode).AND.(iprint_DM>=2)) &
             !        write(io_lun,*) 'Tolerance achieved - exiting'
             !   done = .true.
             !   ndone = n_iter
             !   return
          endif
          call stop_print_timer(tmr_l_iter, "an earlyDM iteration", &
               IPRINT_TIME_THRES1)
       enddo ! End of main loop
    end if ! if (flag_spin_polarisation) then


    if (flag_spin_polarisation) then
       call free_temp_matrix(mat_search_dn)
       call free_temp_matrix(matSphi_dn)
       call free_temp_matrix(matSM3_dn)
       call free_temp_matrix(matM3_dn)
    end if
    call free_temp_matrix(mat_search)
    call free_temp_matrix(mat_temp)
    call free_temp_matrix(matSphi)
    call free_temp_matrix(matSM3)
    call free_temp_matrix(matM3)

1   format('Iteration: ',i3,' Energy: ',e20.12,' Residual: ',e20.12)
2   format('Change in energy: ',e20.12,' Linearity: ',e20.12)

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
!!    Added spin polarisation
!!    Removed redundant parameter number_of_bands
!!  SOURCE
!!
  subroutine lateDM (ndone, n_L_iterations, done, deltaE, vary_mu, mu,&
       inode, ionode, tolerance, record)

    use datatypes
    use numbers
    use logicals
    use matrix_data, ONLY: Lrange, TLrange
    use mult_module, ONLY: matrix_product, allocate_temp_matrix, &
         free_temp_matrix, LNV_matrix_multiply, matphi, matphi_dn, &
         matT, matL, matL_dn, mult, T_L_TL, TL_T_L, matrix_sum, &
         matrix_product_trace, matrix_scale
    use Pulay
    use PosTan, ONLY: PulayR, PulayE, max_iters
    use GenComms, ONLY: cq_abort, gsum
    use global_module, ONLY: iprint_DM,IPRINT_TIME_THRES1,ni_in_cell, &
         flag_global_tolerance, flag_mix_L_SC_min, &
         flag_spin_polarisation, flag_fix_spin_population
    use timer_module, ONLY: cq_timer,start_timer, stop_print_timer, &
         WITH_LEVEL
    use io_module, ONLY: dump_matrix
    use functions_on_grid, ONLY: supportfns, H_on_supportfns
    use H_matrix_module, ONLY: get_H_matrix
    use density_module, ONLY: density, density_up, density_dn, &
         get_electronic_density
    use maxima_module, ONLY: maxngrid
   !Prints out charge density -- 2010.Nov.06 TM
    use io_module, ONLY: dump_charge
    use dimens, ONLY: n_my_grid_points

    implicit none

    ! Passed variables
    integer :: n_L_iterations, inode, ionode, length, ndone
    logical :: vary_mu, done, record
    real(double) :: mu, tolerance, deltaE

    ! Local variables
    integer, dimension(maxpulayDMM) :: mat_Lstore, mat_Gstore,&
         & mat_SGstore
    integer, dimension(maxpulayDMM) :: mat_Lstore_dn, mat_Gstore_dn,&
         & mat_SGstore_dn
    integer :: matM3, matSM3, matSphi, mat_temp
    integer :: matM3_dn, matSM3_dn, matSphi_dn
    real(double) :: e_dot_n, e_dot_n_up, e_dot_n_dn
    real(double) :: n_dot_n, n_dot_n_up, n_dot_n_dn
    real(double) :: electrons, electrons_up, electrons_dn
    real(double) :: energy0, energy0_up, energy0_dn
    real(double) :: energy1, energy1_up, energy1_dn
    real(double) :: g0, g1, gg, step, dsum_factor
    real(double) :: Aij(maxpulayDMM, maxpulayDMM), alph(maxpulayDMM)
    real(double) :: Aij1(maxpulayDMM*maxpulayDMM)
    integer :: n_iter, i,j, pul_mx, npmod
    type(cq_timer) :: tmr_l_tmp1,tmr_l_iter
    !TM
    integer :: iter_stuck = 0
    integer,parameter :: mx_stuck = 5

    iter_stuck=0

    call start_timer(tmr_l_tmp1,WITH_LEVEL)
    if (ndone > n_L_iterations) &
         call cq_abort('lateDM: too many L iterations', &
                       ndone, n_L_iterations)
    do i = 1, maxpulayDMM
       mat_Lstore(i) = allocate_temp_matrix(Lrange,0)
       mat_Gstore(i) = allocate_temp_matrix(Lrange,0)
       mat_SGstore(i) = allocate_temp_matrix(Lrange,0)
    end do
    matM3 = allocate_temp_matrix(Lrange,0)
    matSM3 = allocate_temp_matrix(Lrange,0)
    matSphi = allocate_temp_matrix(Lrange,0)
    mat_temp = allocate_temp_matrix(TLrange,0)

    if (flag_spin_polarisation) then

       ! allocate the required extra matrices and arrays
       do i=1,maxpulayDMM
          mat_Lstore_dn(i) = allocate_temp_matrix(Lrange,0)
          mat_Gstore_dn(i) = allocate_temp_matrix(Lrange,0)
          mat_SGstore_dn(i) = allocate_temp_matrix(Lrange,0)
       end do
       matM3_dn = allocate_temp_matrix(Lrange,0)
       matSM3_dn = allocate_temp_matrix(Lrange,0)
       matSphi_dn = allocate_temp_matrix(Lrange,0)
       ! Update the charge density if flag is set
       if (flag_mix_L_SC_min) then
          ! 2011/08/29 L.Tong: 
          ! original also calculates matphi, (dophi), but I think
          ! this is redundant. So set dontphi. Only doK.
          call LNV_matrix_multiply (electrons_up, energy1_up, doK, &
               dontM1, dontM2, dontM3, dontM4, dontphi, doE, 0, 0, 0, &
               0, spin=1)
          call LNV_matrix_multiply (electrons_dn, energy1_dn, doK, &
               dontM1, dontM2, dontM3, dontM4, dontphi, doE, 0, 0, 0, 0,&
               spin=2)
          energy1 = energy1_up + energy1_dn
          ! note H_on_supportfns is used just as a temp working array
          call get_electronic_density (density_up, electrons_up, &
               supportfns, H_on_supportfns, inode, ionode, maxngrid, &
               spin=1)
          call get_electronic_density (density_dn, electrons_dn, &
               supportfns, H_on_supportfns, inode, ionode, maxngrid, &
               spin=2)
          call get_H_matrix (.true., .false., electrons_up, &
               electrons_dn, density_up, density_dn, maxngrid)
       end if
       ! Get the gradient at the starting point (?) this updates matM3
       ! and matphi
       call LNV_matrix_multiply (electrons_up, energy0_up, dontK, &
            dontM1, dontM2, doM3, dontM4, dophi, doE, 0, matM3, 0, &
            matphi, spin=1)
       call LNV_matrix_multiply (electrons_dn, energy0_dn, dontK, &
            dontM1, dontM2, doM3, dontM4, dophi, doE, 0, matM3_dn, 0, &
            matphi_dn, spin= 2)
       electrons = electrons_up + electrons_dn
       energy0 = energy0_up + energy0_dn
       ! Covariant gradient in SM3
       call matrix_product (matT, matM3, mat_temp, mult(T_L_TL))
       call matrix_product (mat_temp, matT, matSM3, mult(TL_T_L))
       call matrix_product (matT, matM3_dn, mat_temp, mult(T_L_TL))
       call matrix_product (mat_temp, matT, matSM3_dn, mult(TL_T_L))
       ! Project electron gradient out
       if (vary_mu) then
          ! update matSphi from matphi
          call matrix_product (matT, matphi, mat_temp, mult(T_L_TL))
          call matrix_product (mat_temp, matT, matSphi, mult(TL_T_L))
          call matrix_product (matT, matphi_dn, mat_temp, mult(T_L_TL))
          call matrix_product (mat_temp, matT, matSphi_dn, mult(TL_T_L))
          e_dot_n_up = matrix_product_trace (matSM3, matphi)
          e_dot_n_dn = matrix_product_trace (matSM3_dn, matphi_dn)
          e_dot_n = e_dot_n_up + e_dot_n_dn
          n_dot_n_up = matrix_product_trace (matSphi, matphi)
          n_dot_n_dn = matrix_product_trace (matSphi_dn, matphi_dn)
          n_dot_n = n_dot_n_up + n_dot_n_dn
          if (inode == ionode .and. iprint_DM >= 2) then
             write (io_lun, *) 'edotn_up, edotn_dn and edotn are: ',&
                               e_dot_n_up, e_dot_n_dn, e_dot_n
             write (io_lun, *) 'ndotn_up, ndotn_dn and ndotn are: ',&
                               n_dot_n_up, n_dot_n_dn, n_dot_n
          end if
          ! matSM3 mat matM3 are used as search direction (-SG and -G)
          ! These should be positive because we SUBTRACT M3 off below
          ! Here, we can alter M3 directly because it's not expected
          ! to be an exact gradient
          if (flag_fix_spin_population) then
             dsum_factor = e_dot_n_up / n_dot_n_up
             call matrix_sum(one, matSM3, -dsum_factor, matSphi)
             call matrix_sum(one, matM3, -dsum_factor, matphi)
             dsum_factor = e_dot_n_dn / n_dot_n_dn
             call matrix_sum(one, matSM3_dn, -dsum_factor, matSphi_dn)
             call matrix_sum(one, matM3_dn, -dsum_factor, matphi_dn)
          else
             dsum_factor = e_dot_n / n_dot_n
             call matrix_sum(one, matSM3, -dsum_factor, matSphi)
             call matrix_sum(one, matM3, -dsum_factor, matphi)
             call matrix_sum(one, matSM3_dn, -dsum_factor, matSphi_dn)
             call matrix_sum(one, matM3_dn, -dsum_factor, matphi_dn)
          end if
       end if
       ! Store initial gradient (these are search directions G and L)
       call matrix_sum (zero, mat_SGstore(1), -one, matSM3)
       call matrix_sum (zero, mat_SGstore_dn(1), -one, matSM3_dn)
       call matrix_sum (zero, mat_Gstore(1), -one, matM3)
       call matrix_sum (zero, mat_Gstore_dn(1), -one, matM3_dn)
       call matrix_sum (zero, mat_Lstore(1), one, matL)
       call matrix_sum (zero, mat_Lstore_dn(1), one, matL_Dn)

    else ! if (flag_spin_polarisation) then

       ! Update the charge density if flag is set
       if (flag_mix_L_SC_min) then
          call LNV_matrix_multiply(electrons, energy1, doK, dontM1, &
               dontM2, dontM3, dontM4, dophi, doE, 0, matM3, 0, &
               matphi)
          call get_electronic_density(density, electrons, supportfns, &
               H_on_supportfns, inode, ionode, maxngrid)
          call get_H_matrix(.true., .false., electrons, density, &
               maxngrid)
       end if
       ! Get the gradient at the starting point (?)
       call LNV_matrix_multiply(electrons, energy0, &
            dontK, dontM1, dontM2, doM3, dontM4, dophi, doE, 0, matM3,&
            0, matphi)
       ! Covariant gradient in SM3
       call matrix_product(matT, matM3, mat_temp, mult(T_L_TL))
       call matrix_product(mat_temp, matT, matSM3, mult(TL_T_L))
       ! Project electron gradient out
       if (vary_mu) then
          call matrix_product(matT, matphi, mat_temp, mult(T_L_TL))
          call matrix_product(mat_temp, matT, matSphi, mult(TL_T_L))
          e_dot_n = matrix_product_trace(matSM3, matphi) 
          n_dot_n = matrix_product_trace(matSphi, matphi)
          if (inode == ionode .and. iprint_DM >= 2) &
               write (io_lun, *) 'edotn, ndotn are: ', &
               e_dot_n, n_dot_n
          ! These should be positive because we SUBTRACT M3 off below
          ! Here, we can alter M3 directly because it's not expected
          ! to be an exact gradient
          call matrix_sum(one, matSM3, -(e_dot_n)/(n_dot_n), matSphi)
          call matrix_sum(one, matM3, -(e_dot_n)/(n_dot_n), matphi)
       end if
       ! Store initial gradient
       call matrix_sum(zero, mat_SGstore(1), -one, matSM3)
       call matrix_sum(zero, mat_Gstore(1), -one, matM3)
       call matrix_sum(zero, mat_Lstore(1), one, matL)

    end if ! if (flag_spin_polarisation) then
    
    ! timer
    call stop_print_timer (tmr_l_tmp1, "lateDM - Preliminaries",&
         & IPRINT_TIME_THRES1)

    !-----------!
    ! MAIN LOOP !
    !-----------!
    if (flag_spin_polarisation) then
       
       if (flag_global_tolerance) then
          g0 = matrix_product_trace(matM3, matSM3)
          g0 = g0 + matrix_product_trace(matM3_dn, matSM3_dn)
       else
          g0 = matrix_product_trace(matM3, matSM3) / &
               real(ni_in_cell, double)
          g0 = g0 + matrix_product_trace(matM3_dn, matSM3_dn) / &
               real(ni_in_cell, double)
       end if
       do n_iter = 1, n_L_iterations
          call start_timer(tmr_l_iter, WITH_LEVEL)
          if (inode == ionode .and. iprint_DM >= 1) &
               write (io_lun, 1) n_iter, energy0, g0
          ! Storage for pulay DMs/residuals
          npmod = mod(n_iter, maxpulayDMM) + 1
          pul_mx = min(n_iter + 1, maxpulayDMM)
          ! Take a step - maybe correct electron number after
          if (flag_global_tolerance) then
             ! Base step on present gradient and expected dE
             step = deltaE / g0 
          else
             ! Base step on present gradient and expected dE
             step = deltaE / (real(ni_in_cell, double) * g0) 
          end if
          ! We don't want the step to be too small or too big
          if (abs(step) < 0.001_double) step = 0.001_double
          if (abs(step) > 0.1_double) step = 0.1_double
          if (inode == ionode .and. iprint_DM >= 2) &
               write (io_lun,*) 'npmod, pul_mx and step: ', &
                                npmod, pul_mx, step
          ! take L_n+1 = L_n + step * G_n
          if (npmod > 1) then
             call matrix_sum (zero, matL, one, mat_Lstore(npmod-1))
             call matrix_sum (zero, matL_dn, one, mat_Lstore_dn(npmod-1))
             call matrix_sum (one, matL, step, mat_SGstore(npmod-1))
             call matrix_sum (one, matL_dn, step, mat_SGstore_dn(npmod-1))
          else
             call matrix_sum (zero, matL, one, mat_Lstore(pul_mx))
             call matrix_sum (zero, matL_dn, one, mat_Lstore_dn(pul_mx))
             call matrix_sum (one, matL, step, mat_SGstore(pul_mx))
             call matrix_sum (one, matL_dn, step, mat_SGstore_dn(pul_mx))
          endif
          ! after the step, correct the electron number
          if (vary_mu) then
             call correct_electron_number (iprint_DM, inode, ionode)
          endif
          ! Re-evaluate gradient and energy
          call LNV_matrix_multiply (electrons_up, energy1_up, dontK, &
               dontM1, dontM2, doM3, dontM4, dophi, doE, 0, matM3, 0, &
               matphi, spin=1)
          call LNV_matrix_multiply (electrons_dn, energy1_dn, dontK, &
               dontM1, dontM2, doM3, dontM4, dophi, doE, 0, matM3_dn, &
               0, matphi_dn, spin=2)
          electrons = electrons_up + electrons_dn
          energy1 = energy1_up + energy1_dn
          ! Covariant gradient in SM3
          call matrix_product (matT, matM3, mat_temp, mult(T_L_TL))
          call matrix_product (mat_temp, matT, matSM3, mult(TL_T_L))
          call matrix_product (matT, matM3_dn, mat_temp, mult(T_L_TL))
          call matrix_product (mat_temp, matT, matSM3_dn, mult(TL_T_L))
          ! Project out electron variation
          if (vary_mu) then
             ! update matSphi
             call matrix_product (matT, matphi, mat_temp, mult(T_L_TL))
             call matrix_product (mat_temp, matT, matSphi, mult(TL_T_L))
             call matrix_product (matT, matphi_dn, mat_temp, mult(T_L_TL))
             call matrix_product (mat_temp, matT, matSphi_dn, mult(TL_T_L))
             e_dot_n_up = matrix_product_trace (matSM3, matphi)
             e_dot_n_dn = matrix_product_trace (matSM3_dn, matphi_dn)
             e_dot_n = e_dot_n_up + e_dot_n_dn
             n_dot_n_up = matrix_product_trace (matSphi, matphi)
             n_dot_n_dn = matrix_product_trace (matSphi_dn, matphi_dn)
             n_dot_n = n_dot_n_up + n_dot_n_dn
             if (flag_fix_spin_population) then
                dsum_factor = e_dot_n_up / n_dot_n_up
                call matrix_sum (one, matSM3, -dsum_factor, matSphi)
                call matrix_sum (one, matM3, -dsum_factor, matphi)
                dsum_factor = e_dot_n_dn / n_dot_n_dn
                call matrix_sum (one, matSM3_dn, -dsum_factor, matSphi_dn)
                call matrix_sum (one, matM3_dn, -dsum_factor, matphi_dn)
             else
                dsum_factor = e_dot_n / n_dot_n
                call matrix_sum (one, matSM3, -dsum_factor, matSphi)
                call matrix_sum (one, matM3, -dsum_factor, matphi)
                call matrix_sum (one, matSM3_dn, -dsum_factor, matSphi_dn)
                call matrix_sum (one, matM3_dn, -dsum_factor, matphi_dn)
             end if
          end if
          ! Find the residual (i.e. the gradient)
          if (flag_global_tolerance) then
             gg = matrix_product_trace (matM3, matSM3)
             gg = gg + matrix_product_trace (matM3_dn, matSM3_dn)
          else
             gg = matrix_product_trace (matM3, matSM3) / &
                  real(ni_in_cell,double)
             gg = gg + matrix_product_trace (matM3_dn, matSM3_dn) / &
                  real(ni_in_cell,double)
          end if
          if (inode == ionode .and. iprint_DM >= 2) &
               write (io_lun, *) 'R2 is ', sqrt(gg)
          ! recored the Pulay histories
          call matrix_sum (zero, mat_SGstore(npmod), -one, matSM3)
          call matrix_sum (zero, mat_SGstore_dn(npmod), -one, matSM3_dn)
          call matrix_sum (zero, mat_Gstore(npmod), -one, matM3)
          call matrix_sum (zero, mat_Gstore_dn(npmod), -one, matM3_dn)
          call matrix_sum (zero, mat_Lstore(npmod), one, matL)
          call matrix_sum (zero, mat_Lstore_dn(npmod), one, matL_dn)
          ! calculate A_ij
          Aij = zero
          do i = 1, pul_mx
             do j = 1, pul_mx
                gg = matrix_product_trace (mat_Gstore(j), mat_SGstore(i))
                gg = gg + &
                     matrix_product_trace (mat_Gstore_dn(j), mat_SGstore_dn(i))
                Aij (j,i) = gg
                Aij1(j + (i-1) * pul_mx) = gg
             enddo
          enddo
          ! Solve to get alphas
          call DoPulay(Aij, alph, pul_mx, maxpulayDMM, inode, ionode)
          ! Make new L matrix from Pulay sum
          call matrix_scale (zero, matL)
          call matrix_scale (zero, matL_dn)
          do i = 1, pul_mx
             call matrix_sum (one, matL, alph(i), mat_Lstore(i))
             call matrix_sum (one, matL_dn, alph(i), mat_Lstore_dn(i))
          enddo
          ! after the step, correct the electron number
          if (vary_mu) then
             call correct_electron_number (iprint_DM, inode, ionode)
          endif
          if (flag_mix_L_SC_min) then
             ! 2011/08/29 L.Tong
             ! the original also has dophi, but I think is is redundant
             call LNV_matrix_multiply (electrons_up, energy1_up, doK, &
                  dontM1, dontM2, dontM3, dontM4, dontphi, doE, 0, 0, &
                  0, 0, spin=1)
             call LNV_matrix_multiply (electrons_dn, energy1_dn, doK, &
                  dontM1, dontM2, dontM3, dontM4, dontphi, doE, 0, 0, &
                  0, 0, spin=2)
             ! H_on_supportfns is used just as a temp working array
             call get_electronic_density (density_up, electrons_up, &
                  supportfns, H_on_supportfns, inode, ionode, &
                  maxngrid, spin=1)
             call get_electronic_density (density_dn, electrons_dn, &
                  supportfns, H_on_supportfns, inode, ionode, &
                  maxngrid, spin=2)
             call get_H_matrix (.true., .false., electrons_up, &
                  electrons_dn, density_up, density_dn, maxngrid)
          end if
          ! re-evaluate the gradient and energy at new position
          call LNV_matrix_multiply (electrons_up, energy1_up, dontK, &
               dontM1, dontM2, doM3, dontM4, dophi, doE, 0, matM3, 0, &
               matphi, spin=1)
          call LNV_matrix_multiply (electrons_dn, energy1_dn, dontK, &
               dontM1, dontM2, doM3, dontM4, dophi, doE, 0, matM3_dn, &
               0, matphi_dn, spin=2)
          energy1 = energy1_up + energy1_dn
          electrons = electrons_up + electrons_dn
          call matrix_product (matT, matM3, mat_temp, mult(T_L_TL))
          call matrix_product (mat_temp, matT, matSM3, mult(TL_T_L))
          call matrix_product (matT, matM3_dn, mat_temp, mult(T_L_TL))
          call matrix_product (mat_temp, matT, matSM3_dn, mult(TL_T_L))
          if (flag_global_tolerance) then
             g1 = matrix_product_trace (matM3, matSM3)
             g1 = g1 + matrix_product_trace (matM3_dn, matSM3_dn)
          else
             g1 = matrix_product_trace (matM3, matSM3) / &
                  real(ni_in_cell, double)
             g1 = g1 + matrix_product_trace (matM3_dn, matSM3_dn) / &
                  real(ni_in_cell, double)
          end if
          if (inode == ionode .and. iprint_DM >= 3) write (io_lun, *) &
               'Residual before electron gradient correction: ', g1
          if (vary_mu) then
             call matrix_product (matT, matphi, mat_temp, mult(T_L_TL))
             call matrix_product (mat_temp, matT, matSphi, mult(TL_T_L))
             call matrix_product (matT, matphi_dn, mat_temp, mult(T_L_TL))
             call matrix_product (mat_temp, matT, matSphi_dn, mult(TL_T_L))
             e_dot_n_up = matrix_product_trace (matSM3, matphi)
             e_dot_n_dn = matrix_product_trace (matSM3_dn, matphi_dn)
             e_dot_n = e_dot_n_up + e_dot_n_dn
             n_dot_n_up = matrix_product_trace (matSphi, matphi)
             n_dot_n_dn = matrix_product_trace (matSphi_dn, matphi_dn)
             n_dot_n = n_dot_n_up + n_dot_n_dn
             if (flag_fix_spin_population) then
                dsum_factor = e_dot_n_up / n_dot_n_up
                call matrix_sum (one, matSM3, -dsum_factor, matSphi)
                call matrix_sum (one, matM3, -dsum_factor, matphi)
                dsum_factor = e_dot_n_dn / n_dot_n_dn
                call matrix_sum (one, matSM3_dn, -dsum_factor, matSphi_dn)
                call matrix_sum (one, matM3_dn, -dsum_factor, matphi_dn)
             else
                dsum_factor = e_dot_n / n_dot_n
                call matrix_sum (one, matSM3, -dsum_factor, matSphi)
                call matrix_sum (one, matM3, -dsum_factor, matphi)
                call matrix_sum (one, matSM3_dn, -dsum_factor, matSphi_dn)
                call matrix_sum (one, matM3_dn, -dsum_factor, matphi_dn)
             end if
          end if
          ! get deltaE
          deltaE = energy1 - energy0
          ! Find the residual
          if (flag_global_tolerance) then
             g1 = matrix_product_trace (matM3, matSM3)
             g1 = g1 + matrix_product_trace (matM3_dn, matSM3_dn)
          else
             g1 = matrix_product_trace (matM3, matSM3)&
                  & / real (ni_in_cell, double)
             g1 = g1 + matrix_product_trace (matM3_dn, matSM3_dn)&
                  & / real (ni_in_cell, double)
          end if
          if (inode == ionode .and. iprint_DM >= 2)&
               & write (io_lun, *) 'New residual: ', g1
          if ((ndone + n_iter) < max_iters) then
             PulayR (ndone + n_iter) = g1
             PulayE (ndone + n_iter) = energy1
          end if
          ! dump the L matrix if required
          if (mod (n_iter, n_dumpL) == 0) then
             call dump_matrix ("L_up", matL, inode)
             call dump_matrix ("L_dn", matL_dn, inode)
          end if
          ! check if tolerance is reached
          if (g1 < tolerance) then 
             done = .true.
             ndone = n_iter
             if (inode == ionode) write (io_lun, *) 'Achieved tolerance in lateDM'
             if (inode == ionode) write (io_lun, fmt='("Final energy&
                  & and residual: ", f24.9, f15.9)') energy1, g1
             call dealloc_lateDM
             call stop_print_timer (tmr_l_iter, "a lateDM iteration",&
                  & IPRINT_TIME_THRES1)
             return
          else if ((.not. flag_mix_L_SC_min) .and. (g1 > two * g0)) then
             if (inode==ionode) write (io_lun, *)&
                  & 'Panic ! Residual increase in lateDM'
             if (inode==ionode) write(io_lun,*)&
                  & 'Final energy and residual: ', energy1, g1
             ndone = n_iter
             call dealloc_lateDM
             call stop_print_timer (tmr_l_iter, "a lateDM&
                  & iteration", IPRINT_TIME_THRES1)
             return
          else if (g1 > 0.99_double * g0) then
             iter_stuck = iter_stuck+1
             if (iter_stuck > mx_stuck) then
                done = .false.
                ndone = n_iter
                if (inode==ionode) write (io_lun,*)&
                     & 'Fail in reducing Residual in lateDM'
                if (inode==ionode) write (io_lun, fmt = '("     &
                     & energy and residual: ",f24.9,f15.9)')&
                     & energy1, g1
                call dealloc_lateDM
                call stop_print_timer (tmr_l_iter, "a lateDM&
                     & iteration", IPRINT_TIME_THRES1)
                return
             endif ! (iter_stuck > mx_stuck) 
          endif
          ! Replace step with that calculated from real L, and prepare for next step
          call matrix_sum (zero, mat_SGstore(npmod), -one, matSM3)
          call matrix_sum (zero, mat_SGstore_dn(npmod), -one, matSM3_dn)
          call matrix_sum (zero, mat_Gstore(npmod), -one, matM3)
          call matrix_sum (zero, mat_Gstore_dn(npmod), -one, matM3_dn)
          call matrix_sum (zero, mat_Lstore(npmod), one, matL)
          call matrix_sum (zero, mat_Lstore_dn(npmod), one, matL_dn)
          g0 = g1
          energy0 = energy1
          call stop_print_timer (tmr_l_iter, "a lateDM iteration",&
               & IPRINT_TIME_THRES1)
       end do
       
       !Prints out charge density
       if (flag_mix_L_SC_min) then
          call dump_charge (density_up, n_my_grid_points, inode, spin=1)
          call dump_charge (density_dn, n_my_grid_points, inode, spin=2)
       endif  ! (flag_mix_L_SC_min) then
       
    else ! if (flag_spin_polarisation) then
       
       if(flag_global_tolerance) then
          g0 = matrix_product_trace(matM3,matSM3)
       else
          g0 = matrix_product_trace(matM3,matSM3)/real(ni_in_cell,double)
       end if
       do n_iter = 1,n_L_iterations
          call start_timer(tmr_l_iter,WITH_LEVEL)
          if(inode==ionode.AND.iprint_DM>=1) write(io_lun,1) n_iter, energy0, g0
          ! Storage for pulay DMs/residuals
          npmod = mod(n_iter, maxpulayDMM)+1
          pul_mx = min(n_iter+1, maxpulayDMM)
          ! Take a step - maybe correct electron number after
          if(flag_global_tolerance) then
             ! Base step on present gradient and expected dE
             step = deltaE/(g0)
          else
             ! Base step on present gradient and expected dE
             step = deltaE/(real(ni_in_cell,double)*g0)
          end if
          ! We don't want the step to be too small or too big
          if(abs(step)<0.001_double) step = 0.001_double
          if(abs(step)>0.1_double) step = 0.1_double
          if(inode==ionode.and.iprint_DM>=2) &
               write(io_lun,*) 'npmod, pul_mx and step: ',npmod, pul_mx,step
          if(npmod>1) then
             call matrix_sum(zero,matL,one,mat_Lstore(npmod-1))
             call matrix_sum(one,matL,step,mat_SGstore(npmod-1))
          else
             call matrix_sum(zero,matL,one,mat_Lstore(pul_mx))
             call matrix_sum(one,matL,step,mat_SGstore(pul_mx))
          endif
          ! after the step, correct the electron number
          if (vary_mu) then
             call correct_electron_number(iprint_DM, inode, ionode)
          endif
          ! Re-evaluate gradient and energy
          call LNV_matrix_multiply(electrons, energy1, dontK, dontM1, &
               dontM2, doM3, dontM4, dophi, doE, 0, matM3, 0, matphi)
          ! Covariant gradient in SM3
          call matrix_product(matT,matM3,mat_temp,mult(T_L_TL))
          call matrix_product(mat_temp,matT,matSM3,mult(TL_T_L))
          ! Project out electron variation
          if (vary_mu) then
             call matrix_product(matT,matphi,mat_temp,mult(T_L_TL))
             call matrix_product(mat_temp,matT,matSphi,mult(TL_T_L))
             e_dot_n = matrix_product_trace(matSM3,matphi)
             n_dot_n = matrix_product_trace(matSphi,matphi)
             call matrix_sum(one,matSM3,-(e_dot_n)/(n_dot_n),matSphi)
             call matrix_sum(one,matM3,-(e_dot_n)/(n_dot_n),matphi)
          end if
          ! Find the residual (i.e. the gradient)
          if(flag_global_tolerance) then
             gg = matrix_product_trace(matM3,matSM3)
          else
             gg = matrix_product_trace(matM3,matSM3)/real(ni_in_cell,double)
          end if
          if(inode==ionode.AND.iprint_DM>=2) write(io_lun,*) 'R2 is ',sqrt(gg)
          call matrix_sum(zero,mat_SGstore(npmod),-one,matSM3)
          call matrix_sum(zero,mat_Gstore(npmod),-one,matM3)
          call matrix_sum(zero,mat_Lstore(npmod),one,matL)
          Aij = zero
          do i=1,pul_mx
             do j=1,pul_mx
                gg = matrix_product_trace(mat_Gstore(j),mat_SGstore(i))
                Aij(j,i) = gg
                Aij1(j+(i-1)*pul_mx) = gg
             enddo
          enddo
          ! Solve to get alphas
          call DoPulay(Aij,alph,pul_mx,maxpulayDMM,inode,ionode)
          ! Make new L matrix
          call matrix_scale(zero,matL)
          do i=1,pul_mx
             call matrix_sum(one,matL,alph(i),mat_Lstore(i))
          enddo
          ! after the step, correct the electron number
          if (vary_mu) then
             call correct_electron_number( iprint_DM, inode, ionode)
          endif
          if(flag_mix_L_SC_min) then
             call LNV_matrix_multiply(electrons, energy1, doK, dontM1,&
                  dontM2, dontM3, dontM4, dophi, doE, 0, matM3, 0, &
                  matphi)
             call get_electronic_density(density, electrons, &
                  supportfns, H_on_supportfns, inode, ionode, &
                  maxngrid)
             call get_H_matrix(.true.,.false.,electrons,density,maxngrid)
          end if
          ! re-evaluate the gradient and energy at new position
          call LNV_matrix_multiply(electrons, energy1, dontK, dontM1, &
               dontM2, doM3, dontM4, dophi, doE, 0, matM3, 0, matphi)
          call matrix_product(matT,matM3,mat_temp,mult(T_L_TL))
          call matrix_product(mat_temp,matT,matSM3,mult(TL_T_L))
          if(flag_global_tolerance) then
             g1 = matrix_product_trace(matM3,matSM3)
          else
             g1 = matrix_product_trace(matM3,matSM3)/real(ni_in_cell,double)
          end if
          if (inode==ionode.and.iprint_DM>=3) &
               write (io_lun,*) &
               'Residual before electron gradient correction: ', g1
          if (vary_mu) then
             call matrix_product(matT,matphi,mat_temp,mult(T_L_TL))
             call matrix_product(mat_temp,matT,matSphi,mult(TL_T_L))
             e_dot_n = matrix_product_trace(matSM3,matphi)
             n_dot_n = matrix_product_trace(matSphi,matphi)
             call matrix_sum(one,matSM3,-(e_dot_n)/(n_dot_n),matSphi)
             call matrix_sum(one,matM3,-(e_dot_n)/(n_dot_n),matphi)
          end if
          deltaE = energy1 - energy0
          ! Find the residual
          if (flag_global_tolerance) then
             g1 = matrix_product_trace(matM3,matSM3)
          else
             g1 = matrix_product_trace(matM3,matSM3)/real(ni_in_cell,double)
          end if
          if (inode==ionode.and.iprint_DM>=2) &
               write(io_lun,*) 'New residual: ',g1
          if (ndone+n_iter<max_iters) then
             PulayR(ndone+n_iter) = g1
             PulayE(ndone+n_iter) = energy1
          end if

          if (mod(n_iter, n_dumpL) == 0) &
               call dump_matrix("L",matL,inode)

          if (g1<tolerance) then 
             done = .true.
             ndone = n_iter
             if (inode==ionode) &
                  write(io_lun,*) 'Achieved tolerance in lateDM'
             if (inode==ionode) &
                  write (io_lun,fmt='("Final energy and residual: ",f24.9,f15.9)')&
                  energy1,g1
             call dealloc_lateDM
             call stop_print_timer(tmr_l_iter, "a lateDM iteration", &
                  IPRINT_TIME_THRES1)
             return
          else if ((.NOT.flag_mix_L_SC_min).AND.(g1>2.0_double*g0)) then
             if (inode==ionode) &
                  write(io_lun,*) 'Panic ! Residual increase in lateDM'
             if (inode==ionode) &
                  write(io_lun,*) 'Final energy and residual: ',energy1,g1
             ndone = n_iter
             call dealloc_lateDM
             call stop_print_timer(tmr_l_iter, "a lateDM iteration", &
                  IPRINT_TIME_THRES1)
             return
             !2011.11.15 TM
             !OLD else if((.NOT.flag_mix_L_SC_min).AND.(g1>0.99_double*g0)) then
          else if (g1 > 0.99_double * g0) then
             iter_stuck = iter_stuck + 1
             if (iter_stuck > mx_stuck) then
                done = .false.
                ndone = n_iter
                if (inode == ionode) &
                     write(io_lun,*) 'Fail in reducing Residual in lateDM'
                if (inode == ionode) &
                     write(io_lun,fmt='("      energy and residual: ",&
                     &f24.9,f15.9)') energy1, g1
                call dealloc_lateDM
                call stop_print_timer(tmr_l_iter,"a lateDM iteration", &
                     IPRINT_TIME_THRES1)
                return
             endif ! (iter_stuck > mx_stuck) 
             !2011.11.15 TM
          endif
          ! Replace step with real L
          call matrix_sum(zero,mat_SGstore(npmod),-one,matSM3)
          call matrix_sum(zero,mat_Gstore(npmod),-one,matM3)
          call matrix_sum(zero,mat_Lstore(npmod),one,matL)
          g0 = g1
          energy0 = energy1
          call stop_print_timer(tmr_l_iter,"a lateDM iteration",IPRINT_TIME_THRES1)
       enddo

       !Commented out : 2010.12.26 TM
       ! if(g1<tolerance*100.0_double) done = .true.
       !Prints out charge density -- 2010.Nov.06 TM
       if(flag_mix_L_SC_min) then
          call dump_charge(density,n_my_grid_points,inode)
       endif  ! (flag_mix_L_SC_min) then
       
    end if ! if (flag_spin_polarisation) then
    
    call dealloc_lateDM
    
1   format('Iteration: ',i3,' Energy: ',e20.12,' Residual: ',e20.12)

    return

   contains
     
     subroutine dealloc_lateDM
       if (flag_spin_polarisation) then
          call free_temp_matrix (matSphi_dn)
          call free_temp_matrix (matSM3_dn)
          call free_temp_matrix (matM3_dn)
          ! must do this in a reverse order
          do i = maxpulayDMM, 1, -1
             call free_temp_matrix (mat_SGstore_dn(i))
             call free_temp_matrix (mat_Gstore_dn(i))
             call free_temp_matrix (mat_Lstore_dn(i))
          end do
       end if
       call free_temp_matrix(mat_temp)
       call free_temp_matrix(matSphi)
       call free_temp_matrix(matSM3)
       call free_temp_matrix(matM3)
       do i=maxpulayDMM,1,-1
          call free_temp_matrix(mat_SGstore(i))
          call free_temp_matrix(mat_Gstore(i))
          call free_temp_matrix(mat_Lstore(i))
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
!!    Fixed (partly ?) definition of linearity and changed to return zeta in interpG
!!   10:09, 13/02/2006 drb 
!!    Removed all explicit references to data_ variables and rewrote in terms of new 
!!    matrix routines
!!   2011/08/23 L.Tong
!!    Fixed a possible bug, if A=0, then the solution for true step should be
!!    truestep = -C/2B instead of C/2B. Because we are solving the equation
!!    3Ax^2 + 2Bx + C = 0, and if A = 0, x = -C/2B.
!!   2012/03/01 L.Tong
!!    renamed to lineMinL_nospin
!!  SOURCE
!!
  subroutine lineMinL_nospin(output_level, matM3, mat_D, mat_temp, &
       matSM3, energy_in, energy_out, delta_e, inode, ionode, inflex, &
       interpG)

    use datatypes
    use numbers
    use logicals
    use matrix_data, ONLY: Lrange
    use mult_module, ONLY: matrix_sum, matrix_product_trace, allocate_temp_matrix, free_temp_matrix, &
         matrix_product, matT, matL, matL_dn, LNV_matrix_multiply, mult, T_L_TL, TL_T_L, symmetrise_L

    implicit none

    ! Passed variables
    integer :: inode, ionode, output_level, matSM3, mat_temp, matM3, mat_D
    logical :: inflex
    real(double) :: delta_e, energy_in, energy_out, interpG

    ! Local variables
    real(double) :: step, truestep, electrons, A, B, C, D, SQ, g0, g1
    real(double) :: energy_0, energy_1, ig0, ig1, igcross, lamG, igcross2, ig1_step, zeta
    integer i,j,k,length
    integer :: matM3old
    
    matM3old = allocate_temp_matrix(Lrange,0)
    ! mat_D is search direction
    call matrix_sum(zero,matM3old,one,matM3)
    g0 = matrix_product_trace(mat_D,matM3)
    ig0 = matrix_product_trace(mat_D,matM3)

    ! We have a cubic in energy with L_out = L + x.D
    ! energy at x = 0 is energy_in...
    energy_0 = energy_in
    ! The step size should depend on expected energy reduction
    step = delta_e/g0
    if(abs(step)<1.0e-2_double) step = 0.01_double
    if(abs(step)>0.1_double) step = 0.1_double
    if(inode==ionode.and.output_level>=2) write(io_lun,*) 'Step is ',step
    ! now we add a step in the direction D
    call matrix_sum(one,matL,step,mat_D)
    call LNV_matrix_multiply(electrons, energy_1, doK, dontM1, dontM2,&
         doM3, dontM4, dontphi, doE,0,matM3,0,0)
    g1 = matrix_product_trace(mat_D,matM3)
    ! Evaluate old and new gradient cross and new gradient magnitude
    igcross = matrix_product_trace(matSM3,matM3)
    call matrix_product(matT, matM3, mat_temp,mult(T_L_TL))
    call matrix_product(mat_temp, matT, matSM3,mult(TL_T_L))       
    ig1 = matrix_product_trace(matSM3,matM3)
    ig1_step = ig1

    ! the cubic is then given by Ax^3+Bx^2+Cx+D with
    D = energy_0
    C = g0
    !ORI B = (three*energy_1 - step*g1 - two*step*C - three*D)/(step*step)
    B = three*(energy_1-D)/(step*step) - (g1 - two*C)/step
    A = (g1 - two*B*step - C)/(three*step*step)

    ! check for -ve square root      
    SQ = four * B * B - 12.0_double * A * C
    if (SQ.LT.0) then
       ! point of inflexion - problem !
       if (inode .eq. ionode) write (*,*) 'Inflexion  approximation:'
       if(inode==ionode) write(io_lun,*) 'A, B, C, D: ',A, B, C, D
       inflex = .true.
       call free_temp_matrix(matM3old)
       return
    else ! Solve cubic
       if (abs(A).gt.1.0d-14) then
          truestep = (-two*B+SQRT(SQ))/(6.0_double*A)
       else
          if (inode .eq. ionode) write (*,*) 'Local quadratic approximation:'
          ! L.Tong shouldn this be -C/(two*B) instead of the orginal
          ! C/(two * B) ?
          ! truestep = C / (two * B)
          truestep = -C/(two*B)
       end if
    end if
    if (INODE.EQ.IONODE.and.output_level>=2) write(io_lun,*) ' STEP : ',truestep
  !TM 09/09/2003
      if(truestep < 1.e-04_double) truestep = 1.e-04_double
  !TM 09/09/2003
    ! This is unreliable
    !energy_out = A*truestep*truestep*truestep + B*truestep*truestep + &
    !     C*truestep + D
    !delta_e = energy_out - energy_in
    ! update value of matL
    call matrix_sum(one,matL,truestep-step,mat_D)

    ! and make sure the result is symmetric: this avoids the creeping in of
    ! asymmetries due to accumulation of errors
    call symmetrise_L( )
    ! Get K, E and M3 before we return
    call LNV_matrix_multiply(electrons, energy_out, &
         doK, dontM1, dontM2, doM3, dontM4, dontphi, doE,0,matM3,0,0)
    delta_e = energy_out - energy_in
    ! Find interpolated gradient for linearity test
    !    lamG = step/truestep
    if(truestep<step) then
       ! We need to compare interpG with g at truestep
       lamG = truestep/step
       interpG = ig0*(one-lamG)*(one-lamG)+ig1*lamG*lamG+ &
            two*igcross*lamG*(one-lamG)
       zeta = (interpG - ig1)/(ig0-ig1)
    else ! truestep is GREATER than step
       ! We need to compare interpG with g at step
       call matrix_product(matT, matM3, mat_temp,mult(T_L_TL))
       call matrix_product(mat_temp, matT, matSM3,mult(TL_T_L))
       ig1 = matrix_product_trace(matSM3,matM3)
       igcross = matrix_product_trace(matSM3,matM3old)
       lamG = step/truestep
       interpG = ig0*(one-lamG)*(one-lamG)+ig1*lamG*lamG+ &
            two*igcross*lamG*(one-lamG)
       zeta = (interpG - ig1_step)/(ig0-ig1_step)
    end if
    interpG = zeta
    if(inode.EQ.ionode.and.output_level>=2) write(io_lun,*) 'energy_1 is ',energy_1
    !    if(INODE.EQ.IONODE) write(io_lun,*) 'diff is ',energy_1-energy_out
    call free_temp_matrix(matM3old)
    return
  end subroutine lineMinL_nospin
!!***


!!****f* DMMin/lineMinL_spin 
!! PURPOSE
!!  Perform line minimisation for spin polarised calculations in the
!!              direction given by mat_D_up \oplus mat_D_dn
!!  based on lineMinL_nospin subroutine
!!
!!  Note that we do line minimisation by moving along the direction of
!!  give search direction mat_D = direct sum of mat_D_up and
!!  mat_D_dn. This is the same for both the case if the spin
!!  population (i.e. total magnetic moment) is fixed or not. The
!!  difference in the two cases only apprear in the choice of the
!!  search direction vectors mat_D_up and mat_D_dn.
!!
!!  Because the minimisation is based on the direct sums of the
!!  auxilary matrices, it is the total energy that will be minimised,
!!  not the individual energy components of different spin
!!  channels. And the step size is a scalar, which is same value for
!!  both spin channels. We do it this way so that the minimisation is
!!  strictly in the direction pointed by mat_D.
!!
!!  mat_D_up and mat_D_dn are the search directions for spin up and
!!  down.
!!
!! INPUTS
!! OUTPUT
!! RETURN VALUE
!! AUTHOR
!!   L.Tong
!! CREATION DATE 
!!   2011/08/16
!! MODIFICATION HISTORY
!! SOURCE
!!
  subroutine lineMinL_spin(output_level, matM3_up, matM3_dn, mat_D_up,&
       mat_D_dn, mat_temp, matSM3_up, matSM3_dn, energy_in, &
       energy_out, delta_e, inode, ionode, inflex, interpG)
    
    use datatypes
    use numbers
    use logicals
    use matrix_data, ONLY: Lrange
    use mult_module, ONLY: matrix_sum, matrix_product_trace,&
         & allocate_temp_matrix, free_temp_matrix, matrix_product,&
         & matT, matL, matL_dn, LNV_matrix_multiply, mult, T_L_TL,&
         & TL_T_L, symmetrise_L
    use global_module, ONLY: flag_fix_spin_population 

    implicit none
    
    ! Passed variables
    integer :: inode, ionode, output_level 
    integer :: matM3_up, matM3_dn, matSM3_up, matSM3_dn, mat_D_up,&
         & mat_D_dn, mat_temp
    logical :: inflex
    real(double) :: delta_e, energy_in, energy_out, interpG
    
    ! Local variables
    real(double) :: step, truestep, directsum_step
    real(double) :: electrons, electrons_up, electrons_dn
    real(double) :: A, B, C, D, SQ
    real(double) :: g0, g1, ig0, ig1, ig1_step, igcross, igcross2
    real(double) :: energy_0
    real(double) :: energy_1, energy_1_up, energy_1_dn
    real(double) :: energy_2_up, energy_2_dn
    real(double) :: lamG, zeta
    integer :: matM3old_up, matM3old_dn
    integer i, j, k, length

    matM3old_up = allocate_temp_matrix (Lrange, 0)
    matM3old_dn = allocate_temp_matrix (Lrange, 0)    
    ! copy matM3(L_old)
    call matrix_sum (zero, matM3old_up, one, matM3_up)
    call matrix_sum (zero, matM3old_dn, one, matM3_dn)
    ! calculate g0 = tr(matD * matM3), trace of direct sum is sum of
    ! traces of spin components
    g0 = matrix_product_trace (mat_D_up, matM3_up)
    g0 = g0 + matrix_product_trace (mat_D_dn, matM3_dn)
    ! calculate ig0
    ! L.Tong: WARNING! I just copied the approach from
    ! lineMinL_nospin, but I think this may be a mistake, and should
    ! be: ig0 = tr(matM3 * matM3), not tr(mat_D, matM3)
    ig0 = matrix_product_trace (mat_D_up, matM3_up)
    ig0 = ig0 + matrix_product_trace (mat_D_dn, matM3_dn)
    ! We have a cubic in energy with L_out = L + x.D
    ! energy at x = 0 is energy_in
    energy_0 = energy_in
    ! The intermediate step size should depend on expected energy reduction
    step = delta_e / g0
    if (abs(step) < 1.0e-2_double) step = 0.01_double
    if (abs (step) > 0.1_double) step = 0.1_double
    if (inode == ionode .and. output_level >= 2) write(io_lun,*)&
         & 'Step is ',step
    ! now we take the step in the direction D
    call matrix_sum (one, matL, step, mat_D_up)
    call matrix_sum (one, matL_dn, step, mat_D_dn)
    ! calcuate matM3, energy and numbers of electrons again for L_step
    call LNV_matrix_multiply (electrons_up, energy_1_up, doK, dontM1,&
         & dontM2, doM3, dontM4, dontphi, doE, 0, matM3_up, 0, 0,&
         & spin=1)
    call LNV_matrix_multiply (electrons_dn, energy_1_dn, doK, dontM1,&
         & dontM2, doM3, dontM4, dontphi, doE, 0, matM3_dn, 0, 0,&
         & spin=2)
    energy_1 = energy_1_up + energy_1_dn
    ! g1 = tr(mat_D(L_old) * matM3(L_step))
    g1 = matrix_product_trace (mat_D_up, matM3_up)
    g1 = g1 + matrix_product_trace (mat_D_dn, matM3_dn)
    ! Evaluate old and new gradient cross and new gradient
    ! magnitude. This is for linear interpolation. igcross =
    ! tr(matSM3(L_old) * matM3(L_step))
    igcross = matrix_product_trace (matSM3_up, matM3_up)
    igcross = igcross + matrix_product_trace (matSM3_dn, matM3_dn)
    ! get the matSM3(L_step)
    call matrix_product (matT, matM3_up, mat_temp, mult(T_L_TL))
    call matrix_product (mat_temp, matT, matSM3_up, mult(TL_T_L))
    call matrix_product (matT, matM3_dn, mat_temp, mult(T_L_TL))
    call matrix_product (mat_temp, matT, matSM3_dn, mult(TL_T_L))
    ! get ig1 = tr(matSM3(L_step) * matM3(L_step)),
    ! this is for linear interpolation
    ig1 = matrix_product_trace (matSM3_up, matM3_up)
    ig1 = ig1 + matrix_product_trace (matSM3_dn, matM3_dn)
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
       if (inode == ionode) write (*, *) 'Inflexion approximation:'
       if (inode == ionode) write (io_lun, *) 'A, B, C, D: ', A, B, C, D
       inflex = .true.
       call free_temp_matrix (matM3old_dn)
       call free_temp_matrix (matM3old_up)
       return
    else ! Solve cubic
       if (abs (A) > 1.0d-14) then
          truestep = (-two * B + SQRT (SQ)) / (six * A)
       else 
          if (inode .eq. ionode) write (*, *) 'Local quadratic approximation:'
          ! L.Tong shouldn this be -C/(two*B) instead of the orginal
          ! C/(two * B) ?
          ! truestep = C / (two * B)
          truestep = - C / (two * B)
       end if
    end if ! if (SQ < 0) then
    if (inode == ionode .and. output_level >= 2) &
         & write(io_lun,*) ' Step : ', truestep
    ! update value of matL
    directsum_step = truestep - step
    call matrix_sum(one, matL, directsum_step, mat_D_up)
    call matrix_sum(one, matL_dn, directsum_step, mat_D_dn)
    ! matrix symmetric: this avoids the creeping in of asymmetries due
    ! to accumulation of errors
    call symmetrise_L (spin = 1)
    call symmetrise_L (spin = 2)
    ! Get K(L_truestep), E(L_truestep) and M3(L_truestep) before we return
    call LNV_matrix_multiply (electrons_up, energy_2_up, doK,&
         & dontM1, dontM2, doM3, dontM4, dontphi, doE, 0, matM3_up, 0, 0, spin=1)
    call LNV_matrix_multiply (electrons_dn, energy_2_dn, doK,&
         & dontM1, dontM2, doM3, dontM4, dontphi, doE, 0, matM3_dn, 0, 0, spin=2)
    energy_out = energy_2_up + energy_2_dn
    delta_e = energy_out - energy_in

    !-----------------------------------------------------------------------
    !  finished calculating energy_out, delta_e, matM3_up, matM3_dn, inflex
    !-----------------------------------------------------------------------

    ! Find interpolated gradient for linearity test (find interpG/zeta)
    if (truestep < step) then
       ! We need to compare interpG with g at truestep
       lamG = truestep / step
       interpG = ig0*(one - lamG)*(one - lamG) +&
            & ig1*lamG*lamG + two*igcross*lamG*(one - lamG)
       zeta = (interpG - ig1) / (ig0 - ig1)
    else ! truestep is GREATER than step
       ! We need to compare interpG with g at step
       ! find matSM3(L_truestep)
       call matrix_product (matT, matM3_up, mat_temp, mult(T_L_TL))
       call matrix_product (mat_temp, matT, matSM3_up, mult(TL_T_L))
       call matrix_product (matT, matM3_dn, mat_temp, mult(T_L_TL))
       call matrix_product (mat_temp, matT, matSM3_dn, mult(TL_T_L))
       ! ig1 is now used to store tr(matSM3(L_truestep) * matM3(L_truestep))
       ! ig1_step stores tr(matSM3(L_step) * matM3(L_step))
       ig1 = matrix_product_trace (matSM3_up, matM3_up)
       ig1 = ig1 + matrix_product_trace (matSM3_dn, matM3_dn)
       ! igcross is now tr(matSM3(L_old) * matM3(L_truestep))
       igcross = matrix_product_trace (matSM3_up, matM3old_up)
       igcross = igcross + matrix_product_trace (matSM3_dn, matM3old_dn)
       lamG = step / truestep
       interpG = ig0*(one - lamG)*(one - lamG) +&
            & ig1*lamG*lamG + two*igcross*lamG*(one - lamG)
       zeta = (interpG - ig1_step) / (ig0 - ig1_step)
    end if
    ! returns interpG as zeta
    interpG = zeta

    if(inode == ionode .and. output_level >= 2) write (io_lun, *)&
         & 'energy_1_up, energy_1_dn, energy_1 are: ', energy_1_up,&
         & energy_1_dn, energy_1

    call free_temp_matrix(matM3old_dn)
    call free_temp_matrix(matM3old_up)

    return
  end subroutine lineMinL_spin
!!*****


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
!! SOURCE
!!
  subroutine correct_electron_number (output_level, inode, ionode)

    use mult_module, only: matL, matL_dn, matphi, matphi_dn
    use global_module, only: flag_spin_polarisation, flag_fix_spin_population

    implicit none
    
    ! Passed variables
    integer :: inode, ionode, output_level

    if (flag_spin_polarisation) then
       if (flag_fix_spin_population) then
          call correct_electron_number_fixspin (output_level, inode, &
               ionode, matL, matphi, spin=1)
          call correct_electron_number_fixspin (output_level, inode, &
               ionode, matL_dn, matphi_dn, spin=2)   
       else
          call correct_electron_number_varspin (output_level, inode, &
               ionode)
       end if
    else
       call correct_electron_number_nospin (output_level, inode, &
            ionode)
    end if
    
    return

  end subroutine correct_electron_number
!!*****


!!****f* DMMin/correct_electron_number_nospin *
!!
!!  NAME
!!   correct_electron_number_nospin
!!  USAGE
!!
!!  PURPOSE
!!   Corrects the L matrix so that the electron number is correct
!!  INPUTS
!!
!!  USES
!!
!!  AUTHOR
!!   E.H.Hernandez
!!  CREATION DATE
!!   10/04/95
!!  MODIFICATION HISTORY
!!   8/11/95 by Chris Goringe
!!    The convergence seems to be very sensitive to the
!!    exact electron number. This routine has been
!!    modified to solve the cubic in electron number
!!    exactly.
!!   19/5/97 DRB to implement HeadGordon
!!   23/6/97 DRB to implement HeadGordon
!!    N.B. Ne = 2Tr[K.S] !
!!   22/05/00 by DRB to use new modules/mat mults
!!   30/05/2001 dave
!!    ROBODoc header
!!   08/06/2001 dave
!!    Changed to use gsum from GenComms
!!   20/06/2001 dave
!!    Included in DMMinModule and removed declaration
!!    of SolveCubic as external
!!   14:23, 26/02/2003 drb
!!    Added check on electron number after linear correction
!!   2004/10/28 drb
!!    Changed so that linear correction is only used for ridiculously small errors
!!   10:09, 13/02/2006 drb
!!    Removed all explicit references to data_ variables and rewrote in terms of new
!!    matrix routines
!!   2011/07/27 L.Tong 
!!    removed the parameter for number_of_bands, this is an obsolete
!!    historical artifact and replaced the correct electron with the
!!    ne_in_cell from global module
!!   Thursday, 2011/08/11 L.Tong
!!    Changed subroutine name to correct_electron_number_nospin
!!    this subroutine is to be used for the case with no spin polarisation
!!  SOURCE
!!
  subroutine correct_electron_number_nospin (output_level, inode, ionode)

    use datatypes
    use logicals
    use numbers
    use matrix_data, ONLY: TLrange, Lrange
    use mult_module, ONLY: matT, matTtran, matphi, matL, &
         allocate_temp_matrix, free_temp_matrix, matrix_product_trace, &
         matrix_product, matrix_sum, mult, T_L_TL, TL_T_L, LNV_matrix_multiply, &
         matrix_transpose
    use primary_module, ONLY: bundle
    use GenComms, ONLY: gsum
    use global_module, ONLY: ne_in_cell

    implicit none

    ! Passed variables
    integer :: inode, ionode, output_level

    ! Local variables
    real(double) :: electrons_0, electrons, energy,  &
         step, step1, electrons2, g0, g1, recA, B, C, D, truestep, dne

    integer :: matTL, matphi2, matSphi, matSphi2
    logical :: done
    integer :: iter
    logical :: do_spin_down
    
    done = .false.
    iter = 0
    electrons_0 = ne_in_cell
    matTL = allocate_temp_matrix(TLrange,0)
    matphi2 = allocate_temp_matrix(Lrange,0)
    matSphi = allocate_temp_matrix(Lrange,0)
    matSphi2 = allocate_temp_matrix(Lrange,0)
    ! get electron number and gradient
    call matrix_transpose(matT, matTtran)

    do while(.NOT.done.AND.(iter<20)) ! Was 20 !

       call LNV_matrix_multiply (electrons, energy, dontK, dontM1, dontM2, &
            dontM3, dontM4, dophi, dontE, 0, 0, 0, matphi)
       if (inode .eq. ionode.and.output_level>=2) write(io_lun,1) electrons       
       call matrix_product(matT, matphi, matTL, mult(T_L_TL))
       call matrix_product(matTL, matTtran, matSphi,mult(TL_T_L))
       g0 = matrix_product_trace(matSphi, matphi)

       ! initial guess is linear correction...
       step1 = ( electrons_0 - electrons) / g0
       step = 0.1_double
       !step = step1
       if (inode.eq.ionode.and.output_level>=2) write(io_lun,*) 'g0, step1 are ',g0,step1
       ! if we are within 0.1% of the correct number, linear will do.
       !if (abs((electrons_0 - electrons)/electrons_0)<1e-6_double) then
       if (abs(electrons_0 - electrons)<1e-9_double) then
          call matrix_sum (one, matL, step1, matSphi)
          !call LNV_matrix_multiply(electrons2,energy, &
          !     dontK, dontM1, dontM2, dontM3, dontM4, dophi, dontE,0,0,0,matphi)
          !if (inode .eq. ionode.and.output_level>=2) write(io_lun,2) electrons2
          call free_temp_matrix(matSphi2)
          call free_temp_matrix(matSphi)
          call free_temp_matrix(matphi2)
          call free_temp_matrix(matTL)
          return
       else
          step = step1
          ! Take a step along the electron gradient
          call matrix_sum (one, matL, step1, matSphi)
          ! now get new electron number
          call LNV_matrix_multiply (electrons2, energy, dontK, dontM1, dontM2, &
               dontM3, dontM4, dophi, dontE, 0, 0, 0, matphi2)
          call matrix_product (matT, matphi2, matTL, mult(T_L_TL))
          call matrix_product (matTL, matTtran, matSphi2, mult(TL_T_L))

          g1 = matrix_product_trace(matphi,matSphi2)
          if (inode.eq.ionode.and.output_level>=2) write(io_lun,*) 'g1, elec2 are ',g1,electrons2

          ! get coefficients of polynomial
          D = electrons - electrons_0
          C = g0
          !ORI B = (3.0_double*electrons2 - g1*step - 2.0_double*g0*step - 3.0_double*electrons) / &
          !ORI     (step*step)
          B = 3.0_double*(electrons2-electrons)/(step*step) - (g1 + 2.0_double*g0)/step
          !ORI recA = (step*step*step) / &
          !ORI      (2.0_double*electrons - 2.0_double*electrons2 + g0*step + g1*step)
          recA = (step*step*step) / &
               (2.0_double*(electrons - electrons2) + (g0 + g1)*step)
          B = B * recA
          C = C * recA
          D = D * recA

          truestep = SolveCubic(B,C,D,step,inode,ionode)
          if(inode.eq.ionode.and.output_level>=2) write(io_lun,*) 'Step, truestep', step,truestep

          call matrix_sum(one,matL,truestep-step,matSphi)
          if(truestep==step)then
             if(inode==ionode) write(io_lun,*) 'Still in linear loop'
             ! check that electron number is correct
             call LNV_matrix_multiply(electrons2,energy, &
                  dontK, dontM1, dontM2, dontM3, dontM4, dophi, dontE,0,0,0,matphi)
             !call  main_matrix_multiply( data_H, data_S, data_L, data_T, &
             !     data_M1a, data_M3a,data_phi2, &
             !     electrons2,energy, &
             !     dontK, dontM1, dontM2, dontM3, dontM4, dophi, dontE  )


             if (inode .eq. ionode.and.output_level>=2) write(io_lun,2) electrons2
             dne = abs(electrons2 - electrons_0)
             !ORI if((dne/electrons_0)<1e-6_double) done = .true.
             if((dne)<1e-4_double) done = .true.
             iter = iter+1
          else
             done = .true.
             call LNV_matrix_multiply(electrons2,energy, &
                  dontK, dontM1, dontM2, dontM3, dontM4, dophi, dontE,0,0,0,matphi)
             !call  main_matrix_multiply( data_H, data_S, data_L, data_T, &
             !     data_M1a, data_M3a,data_phi2, &
             !     electrons2,energy, &
             !     dontK, dontM1, dontM2, dontM3, dontM4, dophi, dontE  )


             if (inode .eq. ionode.and.output_level>=2) write(io_lun,2) electrons2
          end if
          ! check that electron number is correct
          !call LNV_matrix_multiply(electrons2,energy, &
          !     dontK, dontM1, dontM2, dontM3, dontM4, dophi, dontE,0,0,0,matphi)
          !if (inode .eq. ionode.and.output_level>=2) write(io_lun,2) electrons2
       endif
    enddo
    call free_temp_matrix(matSphi2)
    call free_temp_matrix(matSphi)
    call free_temp_matrix(matphi2)
    call free_temp_matrix(matTL)

1   format(/20x,'Electron number before correction :',f30.15)
2   format(20x,'Electron number after correction  :',f30.15/)

    return
  end subroutine correct_electron_number_nospin
  !!***


!!****f* DMMin/correct_electron_number_fixspin
!! PURPOSE
!!   correct_electron_numbers for each spin component, treating
!!   the spin populations fixed
!!
!!   It takes mat_phi as an input so that the same subroutine can be
!!   used for both spin components
!! INPUTS
!! OUTPUT
!! RETURN VALUE
!! AUTHOR
!!   L.Tong
!! CREATION DATE 
!!   2011/08/11
!! MODIFICATION HISTORY
!! SOURCE
!!
  subroutine correct_electron_number_fixspin (output_level, inode, &
       ionode, mat_L, mat_phi, spin)
    
    use datatypes
    use logicals
    use numbers
    use matrix_data, ONLY: TLrange, Lrange
    use mult_module, ONLY: matT, matTtran, allocate_temp_matrix, &
         free_temp_matrix, matrix_product_trace, matrix_product, &
         matrix_sum, mult, T_L_TL, TL_T_L, LNV_matrix_multiply, &
         matrix_transpose
    use primary_module, ONLY: bundle
    use GenComms, ONLY: gsum, cq_abort
    use global_module, ONLY: ne_up_in_cell, ne_dn_in_cell

    implicit none

    ! Passed variables
    integer :: inode, ionode, output_level, spin
    integer :: mat_L, mat_phi

    ! Local variables
    real(double) :: electrons_0, electrons, energy, step, step1, &
         electrons_2, g0, g1, recA, B, C, D, truestep, dne
    integer :: matTL, matphi2, matSphi, matSphi2
    logical :: done
    integer :: iter

    done = .false.
    iter = 0
    ! set electrons_0 to the correct electron number
    if (spin == 1) then
       electrons_0 = ne_up_in_cell
    else if (spin == 2) then
       electrons_0 = ne_dn_in_cell
    else
       call cq_abort ('correct_electron_numbers_fixspin: spin must be &
            &1 or 2. spin = ', spin)
    end if
    ! allocate temporary matrices
    matTL = allocate_temp_matrix (TLrange,0)
    matphi2 = allocate_temp_matrix (Lrange,0)
    matSphi = allocate_temp_matrix (Lrange,0)
    matSphi2 = allocate_temp_matrix (Lrange,0)
    ! get electron number and gradient
    call matrix_transpose (matT, matTtran)
    do while ( .not. done .and. (iter < 20))
       call LNV_matrix_multiply (electrons, energy, dontK, dontM1, dontM2, &
            dontM3, dontM4, dophi, dontE, 0, 0, 0, mat_phi, spin)
       if (inode .eq. ionode .and. output_level >= 2) then
          if (spin == 1) then
             write (io_lun,*) 'Spin up electron number before correction', &
                  electrons
          else
             write (io_lun,*) 'Spin down electron number before correction', &
                  electrons
          end if
       end if
       call matrix_product (matT, mat_phi, matTL, mult(T_L_TL))
       call matrix_product (matTL, matTtran, matSphi, mult(TL_T_L))
       g0 = matrix_product_trace (matSphi, mat_phi)
       ! initial guess is linear correction...
       step1 = (electrons_0 - electrons) / g0
       step = 0.1_double
       if (inode == ionode .and. output_level >= 2) &
            write (io_lun, *) 'g0, step1 are ', g0, step1
       if (abs (electrons_0 - electrons) < 1e-9_double) then
          call matrix_sum (one, mat_L, step1, matSphi)
          ! check that electron number is correct
          if (output_level >= 2) then
             call LNV_matrix_multiply (electrons_2, energy, dontK, &
                  dontM1, dontM2, dontM3, dontM4, dophi, dontE, 0, 0, &
                  0, mat_phi, spin)
             if (inode == ionode) then
                if (spin == 1) then
                   write (io_lun, *) &
                        'Spin up electron number after correction', &
                        electrons_2
                else
                   write (io_lun, *) &
                        'Spin down electron number after correction', &
                        electrons_2                   
                end if
             end if
          end if
          call free_temp_matrix (matSphi2)
          call free_temp_matrix (matSphi)
          call free_temp_matrix (matphi2)
          call free_temp_matrix (matTL)
          return
       else
          step = step1
          ! Take a step along the electron gradient
          call matrix_sum (one, mat_L, step1, matSphi)
          call LNV_matrix_multiply (electrons_2, energy, dontK, dontM1, dontM2, &
               dontM3, dontM4, dophi, dontE, 0, 0, 0, matphi2, spin)
          call matrix_product (matT, matphi2, matTL, mult(T_L_TL))
          call matrix_product (matTL, matTtran, matSphi2, mult(TL_T_L))
          g1 = matrix_product_trace (mat_phi, matSphi2)
          if (inode.eq.ionode.and.output_level>=2) &
               write(io_lun,*) 'g1, elec2 are ', g1, electrons_2
          ! get coefficients of polynomial
          D = electrons - electrons_0
          C = g0
          B = three * (electrons_2 - electrons) / (step * step) - &
               (g1 + two * g0) / step
          recA = (step * step * step) / &
               (two * (electrons - electrons_2) + (g0 + g1) * step)
          B = B * recA
          C = C * recA
          D = D * recA
          truestep = SolveCubic (B, C, D, step, inode, ionode)
          if (inode .eq. ionode .and. output_level >= 2) &
               write (io_lun, *) 'Step, truestep', step, truestep
          call matrix_sum (one, mat_L, truestep - step, matSphi)
          if (truestep == step) then
             if (inode == ionode) write (io_lun, *) 'Still in linear loop'
             ! check that electron number is correct
             call LNV_matrix_multiply (electrons_2, energy, dontK, &
                  dontM1, dontM2, dontM3, dontM4, dophi, dontE, 0, 0, &
                  0, mat_phi, spin)
             if (inode == ionode .and. output_level >= 2) then
                if (spin == 1) then
                   write (io_lun, *) &
                        'Spin up electron number after correction', &
                        electrons_2
                else
                   write (io_lun, *) &
                        'Spin down electron number after correction', &
                        electrons_2                   
                end if
             end if
             dne = abs (electrons_2 - electrons_0)
             if ((dne) < 1e-4_double) done = .true.
             iter = iter + 1
          else
             done = .true.
             call LNV_matrix_multiply (electrons_2, energy, dontK, &
                  dontM1, dontM2, dontM3, dontM4, dophi, dontE, 0, 0, &
                  0, mat_phi, spin)
             if (inode == ionode .and. output_level >= 2) then
                if (spin == 1) then
                   write (io_lun, *) &
                        'Spin up electron number after correction', &
                        electrons_2
                else
                   write (io_lun, *) &
                        'Spin down electron number after correction', &
                        electrons_2                   
                end if
             end if
          end if
       end if
    end do
    call free_temp_matrix (matSphi2)
    call free_temp_matrix (matSphi)
    call free_temp_matrix (matphi2)
    call free_temp_matrix (matTL)
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
!!  SOURCE
!!
  subroutine correct_electron_number_varspin (output_level, inode, ionode)

    use datatypes
    use logicals
    use numbers
    use matrix_data, ONLY: TLrange, Lrange
    use mult_module, ONLY: matT, matTtran, matphi, matL, matphi_dn, &
         matL_dn, allocate_temp_matrix, free_temp_matrix, &
         matrix_product_trace, matrix_product, matrix_sum, mult, &
         T_L_TL, TL_T_L, LNV_matrix_multiply, matrix_transpose
    use primary_module, ONLY: bundle
    use GenComms, ONLY: gsum
    use global_module, ONLY: ne_in_cell

    implicit none

    ! Passed variables
    integer :: inode, ionode, output_level

    ! Local variables
    real(double) :: electrons_0, electrons, energy, step, step1, &
         electrons2, g0, g1, recA, B, C, D, truestep, dne
    real(double) :: electrons_up, electrons_down, energy_up, &
         energy_down
    integer :: matTL, matphi2, matSphi, matSphi2
    integer :: matphi2_dn, matSphi_dn, matSphi2_dn
    logical :: done
    integer :: iter

    done = .false.
    iter = 0

    ! set the correct electron number
    !    electrons_0 = two*number_of_bands
    electrons_0 = ne_in_cell

    matTL = allocate_temp_matrix (TLrange,0)
    matphi2 = allocate_temp_matrix (Lrange,0)
    matSphi = allocate_temp_matrix (Lrange,0)
    matSphi2 = allocate_temp_matrix (Lrange,0)
    matphi2_dn = allocate_temp_matrix (Lrange,0)
    matSphi_dn = allocate_temp_matrix (Lrange,0)
    matSphi2_dn = allocate_temp_matrix (Lrange,0)

    ! get electron number and gradient
    call matrix_transpose (matT, matTtran)
    do while (.NOT.done.AND.(iter<20)) ! Was 20 !
       call LNV_matrix_multiply (electrons_up, energy_up, dontK, &
            dontM1, dontM2, dontM3, dontM4, dophi, dontE, 0, 0, 0, &
            matphi, spin=1)
       call LNV_matrix_multiply (electrons_down, energy_down, dontK, &
            dontM1, dontM2, dontM3, dontM4, dophi, dontE, 0, 0, 0, &
            matphi_dn, spin=2)
       electrons = electrons_up + electrons_down
       if (inode .eq. ionode .and. output_level >= 2) &
            write (io_lun, 1) electrons_up, electrons_down, electrons
       call matrix_product (matT, matphi, matTL, mult(T_L_TL))
       call matrix_product (matTL, matTtran, matSphi, mult(TL_T_L))
       g0 = matrix_product_trace (matSphi, matphi)
       ! added second spin contribution if required
       call matrix_product (matT, matphi_dn, matTL, mult(T_L_TL))
       call matrix_product (matTL, matTtran, matSphi_dn, mult(TL_T_L))
       g0 = g0 + matrix_product_trace (matSphi_dn, matphi_dn)

       ! initial guess is linear correction...
       step1 = (electrons_0 - electrons) / g0
       step = 0.1_double
       !step = step1
       if (inode == ionode .and. output_level >= 2) &
            write (io_lun, *) 'g0, step1 are ', g0, step1
       ! if we are within 0.1% of the correct number, linear will do.
       if (abs (electrons_0 - electrons) < 1e-9_double) then
          call matrix_sum (one, matL, step1, matSphi)
          call matrix_sum (one, matL_dn, step1, matSphi_dn)
          ! finish the calculations and free temporary memory
          call free_temp_matrix (matSphi2_dn)
          call free_temp_matrix (matSphi_dn)
          call free_temp_matrix (matphi2_dn)
          call free_temp_matrix (matSphi2)
          call free_temp_matrix (matSphi)
          call free_temp_matrix (matphi2)
          call free_temp_matrix (matTL)
          return
       else
          step = step1
          ! Take a step along the electron gradient
          call matrix_sum (one, matL, step1, matSphi)
          call matrix_sum (one, matL_dn, step1, matSphi_dn)
          call LNV_matrix_multiply (electrons_up, energy_up, dontK, &
               dontM1, dontM2, dontM3, dontM4, dophi, dontE, 0, 0, 0, &
               matphi2, spin=1)
          call LNV_matrix_multiply (electrons_down, energy_down, dontK, &
               dontM1, dontM2, dontM3, dontM4, dophi, dontE, 0, 0, 0, &
               matphi2_dn, spin=2)
          electrons2 = electrons_up + electrons_down
          call matrix_product (matT, matphi2, matTL, mult(T_L_TL))
          call matrix_product (matTL, matTtran, matSphi2, mult(TL_T_L))
          g1 = matrix_product_trace (matphi, matSphi2)
          call matrix_product (matT, matphi2_dn, matTL, mult(T_L_TL))
          call matrix_product (matTL, matTtran, matSphi2_dn, mult(TL_T_L))
          g1 = g1 + matrix_product_trace (matphi_dn, matSphi2_dn)
          if (inode == ionode .and. output_level >= 2) &
               write (io_lun, *) &
               'g1, elec_up, elec_dn, elec2 are', &
               g1, electrons_up, electrons_down, electrons2

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

          truestep = SolveCubic (B, C, D, step, inode, ionode)
          if (inode == ionode .and. output_level >= 2) &
               write (io_lun, *) 'Step, truestep', step, truestep

          call matrix_sum (one, matL, truestep - step, matSphi)
          call matrix_sum (one, matL_dn, truestep - step, matSphi_dn)

          if (truestep == step) then
             if (inode == ionode) write (io_lun, *) 'Still in linear loop'
             ! check that electron number is correct
             call LNV_matrix_multiply (electrons_up, energy_up, dontK,&
                  dontM1, dontM2, dontM3, dontM4, dophi, dontE, 0, 0, &
                  0, matphi, spin=1)
             call LNV_matrix_multiply (electrons_down, energy_down, dontK,&
                  dontM1, dontM2, dontM3, dontM4, dophi, dontE, 0, 0, &
                  0, matphi_dn, spin=2)
             electrons2 = electrons_up + electrons_down
             if (inode == ionode .and. output_level >= 2) &
                  write (io_lun, 2) electrons_up, electrons_down, electrons2
             dne = abs (electrons2 - electrons_0)
             if ((dne / electrons_0) < 1e-6_double) done = .true.
             iter = iter + 1
          else
             done = .true.
             call LNV_matrix_multiply (electrons_up, energy_up, dontK,&
                  dontM1, dontM2, dontM3, dontM4, dophi, dontE, 0, 0, &
                  0, matphi, spin=1)
             call LNV_matrix_multiply (electrons_down, energy_down, dontK,&
                  dontM1, dontM2, dontM3, dontM4, dophi, dontE, 0, 0, &
                  0, matphi_dn, spin=2)
             electrons2 = electrons_up + electrons_down
             if (inode == ionode .and. output_level >= 2) &
                  write (io_lun, 2) electrons_up, electrons_down, electrons2
          end if
       end if
    enddo

    ! free temp matrices
    call free_temp_matrix (matSphi2_dn)
    call free_temp_matrix (matSphi_dn)
    call free_temp_matrix (matphi2_dn)
    call free_temp_matrix (matSphi2)
    call free_temp_matrix (matSphi)
    call free_temp_matrix (matphi2)
    call free_temp_matrix (matTL)

1   format (20x, 'Electron numbers before correction (Nup, Ndn, Ntotal) :', &
         f30.15, f30.15, f30.15)
2   format (20x, 'Electron number after correction (Nup, Ndn, Ntotal) :', &
         f30.15, f30.15, f30.15)

    return
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
!!   Taken from NR, 2nd Edition (p179)
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
!!  SOURCE
!!
  function SolveCubic(A,B,C,guess,inode,ionode)

    use datatypes
    use numbers

    implicit none

    integer inode,ionode
    real(double) :: SolveCubic
    real(double) :: A,B,C,guess

    real(double) :: Q, R, theta, T, Z1, Z2, Z3
    real(double) :: S1, S2

    Q = (A*A-3.0_double*B)/9.0_double
    R = (2.0_double*A*A*A-9.0_double*A*B+27.0_double*C)/54.0_double

    if ((R*R)<(Q*Q*Q)) then
       ! three roots...
       theta = acos(R/(SQRT(Q*Q*Q)))
       T = -2.0_double*SQRT(Q)
       Z1 = T * cos(theta/3.0_double) - (A/3.0_double) - guess
       Z2 = T * cos((theta+2.0_double*pi)/3.0_double) - (A/3.0_double) - guess
       Z3 = T * cos((theta-2.0_double*pi)/3.0_double) - (A/3.0_double) - guess

       if (ABS(Z1)<=ABS(Z2).and.ABS(Z1)<=ABS(Z3)) then
          SolveCubic = Z1 + guess
       else if (ABS(Z2)<=ABS(Z3)) then
          SolveCubic = Z2 + guess
       else
          SolveCubic = Z3 + guess
       end if
    else

       ! the following is based on Abramowitz-Stegun p.17, and it works
       Q = -Q
       R = -R

       S1 = R + sqrt(Q*Q*Q + R*R)
       S2 = R - sqrt(Q*Q*Q + R*R)

       if ( S1 >= zero ) then
          S1 = S1**(1.0_double/3.0_double)
       else
          S1 = -(abs(S1))**(1.0_double/3.0_double)
       end if

       if ( S2 >= zero ) then
          S2 = S2**(1.0_double/3.0_double)
       else
          S2 = -(abs(S2))**(1.0_double/3.0_double)
       end if

       SolveCubic = S1 + S2 - A/3.0_double
       if(abs(SolveCubic)>10.0_double*abs(guess)) then
          if(inode==ionode) write(io_lun,*) 'Step too large: linear guess was: ',SolveCubic, guess
          SolveCubic = guess
       end if

    end if

    return
  end function SolveCubic
!!***
end module DMMin
