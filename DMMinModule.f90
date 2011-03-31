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
!!  SOURCE
!!
  subroutine FindMinDM(n_L_iterations, number_of_bands,  &
       vary_mu, tolerance, mu, inode, ionode, resetL, record)

    use datatypes
    use numbers
    use global_module, ONLY: iprint_DM,IPRINT_TIME_THRES1
    use mult_module, ONLY: matrix_transpose, matT, matTtran, matL, matK, matrix_sum
    use McWeeny, ONLY: InitMcW, McWMin
    use PosTan, ONLY: max_iters, cscale, PulayE, PulayR, &
         PulayC, PulayBeta, pos_tan, fit_coeff
    use DiagModule, ONLY: FindEvals, diagon
    use io_module, ONLY: dump_matrix
    use energy, ONLY: entropy
    use timer_module, ONLY: cq_timer,start_timer,stop_print_timer,WITH_LEVEL

    implicit none

    ! Passed variables
    integer :: n_L_iterations, inode, ionode
    real(double) :: number_of_bands, mu, tolerance
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
       call FindEvals(2.0_double*number_of_bands)
       call stop_timer(tmr_std_densitymat)
       return
    end if
    if(inode==ionode) write(io_lun,*) 'Welcome to FindDMMin, tol: ',tolerance, n_L_iterations
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
          call InitMcW(number_of_bands,inode,ionode)
          call McWMin( n_L_iterations, delta_e ) 
          early = .true.
          resetL = .false.    !2010.Nov.06 TM
       endif
       if(early.OR.problem) then
          inflex = .false.
          call earlyDM(ndone, n_L_iterations, number_of_bands,  &
               delta_e, done, vary_mu, tolerance, mu, &
               inode, ionode, inflex, record)
       endif
       if((.NOT.done).AND.(.NOT.inflex)) then ! Continue on to late stage
          call lateDM(ndone, n_L_iterations, done, delta_e, &
               vary_mu, mu, inode, ionode, number_of_bands, tolerance, record)
       endif
       niter=niter+1
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
!!  SOURCE
!!
  subroutine earlyDM(ndone, n_L_iterations, &
       number_of_bands,delta_e,done,vary_mu, tolerance, mu, &
       inode, ionode, inflex, record)

    use datatypes
    use logicals
    use numbers
    use matrix_data, ONLY: Lrange, TLrange
    use mult_module, ONLY: matT, matphi, mult, T_L_TL, TL_T_L, matrix_product_trace, &
         LNV_matrix_multiply, allocate_temp_matrix, free_temp_matrix, matrix_sum, matrix_product
    use primary_module, ONLY: bundle
    use PosTan, ONLY: PulayR, PulayE
    use GenComms, ONLY: cq_abort, gsum
    use global_module, ONLY: iprint_DM,IPRINT_TIME_THRES1,ni_in_cell,flag_global_tolerance
    use timer_module, ONLY: cq_timer,start_timer, stop_print_timer, WITH_LEVEL

    implicit none

    ! Passed variables
    real(double) ::  number_of_bands, mu, tolerance, delta_e
    integer :: inode, ionode, n_L_iterations, ndone
    logical :: inflex, vary_mu, done, record

    ! Local variables
    integer :: n_iter, length
    real(double) :: g0, interpG, g1, zeta
    real(double) :: e_dot_n, n_dot_n, energy0, energy1, electrons
    integer :: matM3, matSM3, matSphi, mat_temp, mat_search
    type(cq_timer) :: tmr_l_tmp1, tmr_l_iter

    call start_timer(tmr_l_tmp1,WITH_LEVEL)
    matM3 = allocate_temp_matrix(Lrange,0)
    matSM3 = allocate_temp_matrix(Lrange,0)
    matSphi = allocate_temp_matrix(Lrange,0)
    mat_temp = allocate_temp_matrix(TLrange,0)
    mat_search = allocate_temp_matrix(Lrange,0)
    if(ndone>n_L_iterations) call cq_abort('earlyDM: too many L iterations', &
         ndone, n_L_iterations)
    ! Set up the gradient and electron number
    if (vary_mu) then
       call correct_electron_number(iprint_DM,number_of_bands,inode,ionode)
    endif
    call LNV_matrix_multiply(electrons, energy0, dontK, dontM1, dontM2, doM3, dontM4, dophi, doE,0,matM3,0,matphi )
    ! Pre- and post-multiply M3 by S^-1 so that it is contravariant
    call matrix_product(matT, matM3, mat_temp,mult(T_L_TL))
    call matrix_product(mat_temp, matT, matSM3,mult(TL_T_L))
    ! Project gradient perpendicular to electron gradient
    call matrix_sum(zero,mat_search,-one,matSM3)
    if(vary_mu) then
       ! Pre- and post-multiply phi by S^-1 so that it is contravariant
       call matrix_product(matT, matphi, mat_temp,mult(T_L_TL))
       call matrix_product(mat_temp, matT, matSphi,mult(TL_T_L))
       ! Only one is pre- and post-multiplied because (A,B) = Tr(ASBS)
       e_dot_n = matrix_product_trace(matSM3,matphi)
       n_dot_n = matrix_product_trace(matSphi,matphi)
       !mu = e_dot_n/n_dot_n
       if(inode.eq.ionode.and.iprint_DM>=2) &
            write(io_lun,*) 'e.n and n.n: ',e_dot_n, n_dot_n
       ! This is right: search dirn is -SM3
       call matrix_sum(one,mat_search,(e_dot_n)/(n_dot_n),matSphi)
       ! Here, we can't alter M3 directly: lineMinL expects REAL gradient
    endif
    call stop_print_timer(tmr_l_tmp1,"earlyDM - Preliminaries",IPRINT_TIME_THRES1)
    !-----------!
    ! MAIN LOOP !
    !-----------!
    delta_e = zero
    do n_iter = 1,n_L_iterations
       call start_timer(tmr_l_iter,WITH_LEVEL)
       ! Gradient before line min
       if(flag_global_tolerance) then
          g0 = matrix_product_trace(matM3,matSM3)
       else
          g0 = matrix_product_trace(matM3,matSM3)/real(ni_in_cell,double)
       end if
       if(inode==ionode.AND.iprint_DM>=1) write(io_lun,1) n_iter, energy0, g0
1      format('Iteration: ',i3,' Energy: ',e20.12,' Residual: ',e20.12)
       call lineMinL(iprint_DM, matM3, mat_search, mat_temp,matSM3,&
            energy0, energy1, delta_e, inode, ionode, inflex,interpG)       
       delta_e = energy1 - energy0
       if(inflex) then
          if(inode==ionode) write(io_lun,*) 'Panic ! Inflexion point found !'
          ndone = n_iter
          call free_temp_matrix(mat_search)
          call free_temp_matrix(mat_temp)
          call free_temp_matrix(matSphi)
          call free_temp_matrix(matSM3)
          call free_temp_matrix(matM3)
          call stop_print_timer(tmr_l_iter,"an earlyDM iteration",IPRINT_TIME_THRES1)
          return !This is a panic sign !
       endif
       ! Gradient after line min - assumes main_matrix_mult at end of lineMinL
       ! Pre- and post-multiply M3 by S^-1 so that it is contravariant
       call matrix_product(matT, matM3, mat_temp,mult(T_L_TL))
       call matrix_product(mat_temp, matT, matSM3,mult(TL_T_L))       
       if(flag_global_tolerance) then
          g1 = matrix_product_trace(matM3,matSM3)
       else
          g1 = matrix_product_trace(matM3,matSM3)/real(ni_in_cell,double)
       end if
       ! Test for linearity or convergence
       !zeta = (interpG - g1)/(g0)
       ! DRB 2004/09/28 Now zeta returned by lineMinL
       zeta = interpG
       if(inode==ionode.AND.iprint_DM>=2) write(io_lun,2) delta_e, zeta
2      format('Change in energy: ',e20.12,' Linearity: ',e20.12)
       if(inode==ionode.AND.iprint_DM>=2) write(io_lun,*) 'zeta...: ',g0,g1,interpG
       if(vary_mu) then
            call correct_electron_number(iprint_DM,number_of_bands,inode,ionode)
       endif
       call LNV_matrix_multiply(electrons, energy0, dontK, dontM1, dontM2, doM3, dontM4, dophi, doE,0,matM3,0,matphi )
       ! Pre- and post-multiply M3 by S^-1 so that it is contravariant
       call matrix_product(matT, matM3, mat_temp,mult(T_L_TL))
       call matrix_product(mat_temp, matT, matSM3,mult(TL_T_L))
       energy0 = energy1
       call matrix_sum(zero,mat_search,-one,matSM3)
       ! Project search direction perpendicular to electron gradient
       if (vary_mu) then
          ! Wrap electron gradient for tensorial correctness
          call matrix_product(matT, matphi, mat_temp,mult(T_L_TL))
          call matrix_product(mat_temp, matT, matSphi,mult(TL_T_L))
          e_dot_n = matrix_product_trace(matSM3,matphi)
          n_dot_n = matrix_product_trace(matSphi,matphi)
          call matrix_sum(one,mat_search,(e_dot_n)/(n_dot_n),matSphi)
          ! Here, we can't alter M3 directly: lineMinL expects REAL gradient
       end if
       if(flag_global_tolerance) then
          g0 = matrix_product_trace(matM3,matSM3)
       else
          g0 = matrix_product_trace(matM3,matSM3)/real(ni_in_cell,double)
       end if
       PulayR(ndone+n_iter) = g0
       PulayE(ndone+n_iter) = energy0
       if(abs(zeta)<LinTol_DMM) then ! We're linear
          ! Commented out test for convergence to force checking of 
          ! convergence by Pulay - DRB, 22/03/01
          ! Added +1 to n_iter for cosmetic reasons
          if(inode==ionode.AND.iprint_DM>=1) write(io_lun,1) n_iter+1, energy0, g0
          if((inode==ionode).AND.(iprint_DM>=2)) &
               write(io_lun,*) 'Linearity satisfied - calling PulayL'
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
          call stop_print_timer(tmr_l_iter,"an earlyDM iteration",IPRINT_TIME_THRES1)
          return
       !else if(abs(delta_E)<tolerance) then
       !   if(inode==ionode.AND.iprint_DM>=1) write(io_lun,1) n_iter, energy0, g0
       !   if((inode==ionode).AND.(iprint_DM>=2)) &
       !        write(io_lun,*) 'Tolerance achieved - exiting'
       !   done = .true.
       !   ndone = n_iter
       !   return
       endif
       call stop_print_timer(tmr_l_iter,"an earlyDM iteration",IPRINT_TIME_THRES1)
    enddo ! End of main loop
    call free_temp_matrix(mat_search)
    call free_temp_matrix(mat_temp)
    call free_temp_matrix(matSphi)
    call free_temp_matrix(matSM3)
    call free_temp_matrix(matM3)
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
!!    Removed all explicit references to data_ variables and rewrote in terms of new 
!!    matrix routines
!!   2008/03/03 18:32 dave
!!    Removed dsqrt
!!   2009/07/08 16:41 dave
!!    Introduced atom-based tolerance (just divide by ni_in_cell)
!!   2010/03/18 14:08 dave
!!    Added code for mixed L and SCF minimisation
!!  SOURCE
!!
  subroutine lateDM(ndone, n_L_iterations, done, deltaE, &
       vary_mu, mu, inode, ionode, number_of_bands, tolerance, record)

    use datatypes
    use numbers
    use logicals
    use matrix_data, ONLY: Lrange, TLrange
    use mult_module, ONLY: matrix_product, allocate_temp_matrix, free_temp_matrix, LNV_matrix_multiply,&
         matphi, matT, matL, mult, T_L_TL, TL_T_L, matrix_sum, matrix_product_trace, matrix_scale
    use Pulay
    use PosTan, ONLY: PulayR, PulayE, max_iters
    use GenComms, ONLY: cq_abort, gsum
    use global_module, ONLY: iprint_DM,IPRINT_TIME_THRES1,ni_in_cell, flag_global_tolerance, flag_mix_L_SC_min
    use timer_module, ONLY: cq_timer,start_timer, stop_print_timer, WITH_LEVEL
    use io_module, ONLY: dump_matrix
    use functions_on_grid, ONLY: supportfns, H_on_supportfns
    use H_matrix_module, ONLY: get_H_matrix
    use density_module, ONLY: density, get_electronic_density
    use maxima_module, ONLY: maxngrid
   !Prints out charge density -- 2010.Nov.06 TM
    use io_module, ONLY: dump_charge
    use dimens, ONLY: n_my_grid_points

    implicit none

    ! Passed variables
    integer :: n_L_iterations, inode, ionode, length, ndone
    logical :: vary_mu, done, record
    real(double) :: mu, number_of_bands, tolerance, deltaE

    ! Local variables
    integer, dimension(maxpulayDMM) :: mat_Lstore, mat_Gstore, mat_SGstore
    integer :: matM3, matSM3, matSphi,mat_temp
    real(double) :: e_dot_n, n_dot_n, electrons, step, energy0, energy1
    real(double) :: g0, g1, gg
    real(double) :: Aij(maxpulayDMM, maxpulayDMM), alph(maxpulayDMM)
    real(double) :: Aij1(maxpulayDMM*maxpulayDMM)
    integer :: n_iter, i,j, pul_mx, npmod
    type(cq_timer) :: tmr_l_tmp1,tmr_l_iter
!TM
    integer :: iter_stuck = 0
    integer,parameter :: mx_stuck = 5

    iter_stuck=0

    call start_timer(tmr_l_tmp1,WITH_LEVEL)
    if(ndone>n_L_iterations) call cq_abort('lateDM: too many L iterations', &
         ndone, n_L_iterations)
    do i=1,maxpulayDMM
       mat_Lstore(i) = allocate_temp_matrix(Lrange,0)
       mat_Gstore(i) = allocate_temp_matrix(Lrange,0)
       mat_SGstore(i) = allocate_temp_matrix(Lrange,0)
    end do
    matM3 = allocate_temp_matrix(Lrange,0)
    matSM3 = allocate_temp_matrix(Lrange,0)
    matSphi = allocate_temp_matrix(Lrange,0)
    mat_temp = allocate_temp_matrix(TLrange,0)
    !if (vary_mu) then
    !   call correct_electron_number( iprint_DM, number_of_bands,inode, ionode)
    !endif
    ! Update the charge density if flag is set
    if(flag_mix_L_SC_min) then
       call LNV_matrix_multiply(electrons, energy1, &
            doK, dontM1, dontM2, dontM3, dontM4, dophi, doE,0,matM3,0,matphi  )
       call get_electronic_density(density, electrons, supportfns, H_on_supportfns, inode, ionode, maxngrid)
       call get_H_matrix(.true.,.false.,electrons,density,maxngrid)
    end if
    ! Get the gradient at the starting point (?)
    call LNV_matrix_multiply(electrons, energy0, &
         dontK, dontM1, dontM2, doM3, dontM4, dophi, doE,0,matM3,0,matphi  )
    ! Covariant gradient in SM3
    call matrix_product(matT,matM3,mat_temp,mult(T_L_TL))
    call matrix_product(mat_temp,matT,matSM3,mult(TL_T_L))
    ! Project electron gradient out
    if (vary_mu) then
       call matrix_product(matT,matphi,mat_temp,mult(T_L_TL))
       call matrix_product(mat_temp,matT,matSphi,mult(TL_T_L))
       e_dot_n = matrix_product_trace(matSM3,matphi) 
       n_dot_n = matrix_product_trace(matSphi,matphi)
       if(inode==ionode.and.iprint_DM>=2) &
            write(io_lun,*) 'edotn, ndotn are: ',e_dot_n, n_dot_n
       ! These should be positive because we SUBTRACT M3 off below
       ! Here, we can alter M3 directly because it's not expected to be an exact gradient
       call matrix_sum(one,matSM3,-(e_dot_n)/(n_dot_n),matSphi)
       call matrix_sum(one,matM3,-(e_dot_n)/(n_dot_n),matphi)
    end if
    ! Store initial gradient
    call matrix_sum(zero,mat_SGstore(1),-one,matSM3)
    call matrix_sum(zero,mat_Gstore(1),-one,matM3)
    call matrix_sum(zero,mat_Lstore(1),one,matL)
    call stop_print_timer(tmr_l_tmp1,"lateDM - Preliminaries",IPRINT_TIME_THRES1)
    !-----------!
    ! MAIN LOOP !
    !-----------!
    if(flag_global_tolerance) then
       g0 = matrix_product_trace(matM3,matSM3)
    else
       g0 = matrix_product_trace(matM3,matSM3)/real(ni_in_cell,double)
    end if
    do n_iter = 1,n_L_iterations
       call start_timer(tmr_l_iter,WITH_LEVEL)
       if(inode==ionode.AND.iprint_DM>=1) write(io_lun,1) n_iter, energy0, g0
1      format('Iteration: ',i3,' Energy: ',e20.12,' Residual: ',e20.12)
       ! Storage for pulay DMs/residuals
       npmod = mod(n_iter, maxpulayDMM)+1
       pul_mx = min(n_iter+1, maxpulayDMM)
       ! Take a step - maybe correct electron number after
       if(flag_global_tolerance) then
          step = deltaE/(g0) ! Base step on present gradient and expected dE
       else
          step = deltaE/(real(ni_in_cell,double)*g0) ! Base step on present gradient and expected dE
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
          call correct_electron_number( iprint_DM, number_of_bands,inode, ionode)
       endif       
       ! Re-evaluate gradient and energy
       call LNV_matrix_multiply(electrons, energy1, &
            dontK, dontM1, dontM2, doM3, dontM4, dophi, doE,0,matM3,0,matphi  )
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
       call DoPulay2D(Aij,alph,pul_mx,maxpulayDMM,inode,ionode)
       ! Make new L matrix
       call matrix_scale(zero,matL)
       do i=1,pul_mx
          call matrix_sum(one,matL,alph(i),mat_Lstore(i))
       enddo
       ! after the step, correct the electron number
       if (vary_mu) then
          call correct_electron_number( iprint_DM, number_of_bands,inode, ionode)
       endif
       if(flag_mix_L_SC_min) then
          call LNV_matrix_multiply(electrons, energy1, &
               doK, dontM1, dontM2, dontM3, dontM4, dophi, doE,0,matM3,0,matphi  )
          call get_electronic_density(density, electrons, supportfns, H_on_supportfns, inode, ionode, maxngrid)
          call get_H_matrix(.true.,.false.,electrons,density,maxngrid)
       end if
       ! re-evaluate the gradient and energy at new position
       call LNV_matrix_multiply(electrons, energy1, &
            dontK, dontM1, dontM2, doM3, dontM4, dophi, doE,0,matM3,0,matphi  )
       call matrix_product(matT,matM3,mat_temp,mult(T_L_TL))
       call matrix_product(mat_temp,matT,matSM3,mult(TL_T_L))
       if(flag_global_tolerance) then
          g1 = matrix_product_trace(matM3,matSM3)
       else
          g1 = matrix_product_trace(matM3,matSM3)/real(ni_in_cell,double)
       end if
       if(inode==ionode.and.iprint_DM>=3) write(io_lun,*) 'Residual before electron gradient correction: ',g1
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
       if(flag_global_tolerance) then
          g1 = matrix_product_trace(matM3,matSM3)
       else
          g1 = matrix_product_trace(matM3,matSM3)/real(ni_in_cell,double)
       end if
       if(inode==ionode.and.iprint_DM>=2) write(io_lun,*) 'New residual: ',g1
       if(ndone+n_iter<max_iters) then
          PulayR(ndone+n_iter) = g1
          PulayE(ndone+n_iter) = energy1
       end if
       
       if(mod(n_iter, n_dumpL) == 0) call dump_matrix("L",matL,inode)

       if(g1<tolerance) then 
          done = .true.
          ndone = n_iter
          if(inode==ionode) write(io_lun,*) 'Achieved tolerance in lateDM'
          if(inode==ionode) write(io_lun,fmt='("Final energy and residual: ",f24.9,f15.9)') energy1,g1
          call dealloc_lateDM
          call stop_print_timer(tmr_l_iter,"a lateDM iteration",IPRINT_TIME_THRES1)
          return
       else if((.NOT.flag_mix_L_SC_min).AND.(g1>2.0_double*g0)) then
          if(inode==ionode) write(io_lun,*) 'Panic ! Residual increase in lateDM'
          if(inode==ionode) write(io_lun,*) 'Final energy and residual: ',energy1,g1
          ndone = n_iter
          call dealloc_lateDM
          call stop_print_timer(tmr_l_iter,"a lateDM iteration",IPRINT_TIME_THRES1)
          return
!2011.11.15 TM
       !OLD else if((.NOT.flag_mix_L_SC_min).AND.(g1>0.99_double*g0)) then
       else if(g1>0.99_double*g0) then
          iter_stuck=iter_stuck+1
         if(iter_stuck > mx_stuck) then
          done = .false.
          ndone = n_iter
          if(inode==ionode) write(io_lun,*) 'Fail in reducing Residual in lateDM'
          if(inode==ionode) write(io_lun,fmt='("      energy and residual: ",f24.9,f15.9)') energy1,g1
          call dealloc_lateDM
          call stop_print_timer(tmr_l_iter,"a lateDM iteration",IPRINT_TIME_THRES1)
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
    call dealloc_lateDM
    return
   contains
    subroutine dealloc_lateDM
     call free_temp_matrix(mat_temp)
     call free_temp_matrix(matSphi)
     call free_temp_matrix(matSM3)
     call free_temp_matrix(matM3)
     do i=maxpulayDMM,1,-1
        call free_temp_matrix(mat_SGstore(i))
        call free_temp_matrix(mat_Gstore(i))
        call free_temp_matrix(mat_Lstore(i))
     end do
    end subroutine dealloc_lateDM
  end subroutine lateDM
!!***

! -----------------------------------------------------------
! Subroutine lineMinL
! -----------------------------------------------------------

!!****f* DMMin/lineMinL *
!!
!!  NAME 
!!   lineMinL
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
!!  SOURCE
!!
  subroutine lineMinL( output_level, matM3, mat_D, mat_temp,matSM3,energy_in, energy_out, delta_e, &
       inode, ionode, inflex, interpG )

    use datatypes
    use numbers
    use logicals
    use matrix_data, ONLY: Lrange
    use mult_module, ONLY: matrix_sum, matrix_product_trace, allocate_temp_matrix, free_temp_matrix, &
         matrix_product, matL, matT, LNV_matrix_multiply, mult, T_L_TL, TL_T_L, symmetrise_L

    implicit none

    ! Passed variables
    integer :: inode, ionode, output_level, matSM3,mat_temp,matM3, mat_D
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
    call LNV_matrix_multiply(electrons, energy_1, doK, dontM1, dontM2, doM3, dontM4, dontphi, doE,0,matM3,0,0)
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
          truestep = C/(two*B)
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
    ! update value of mat_L
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
  end subroutine lineMinL
!!***

! -----------------------------------------------------------
! Subroutine correct_electron_number
! -----------------------------------------------------------

!!****f* DMMin/correct_electron_number *
!!
!!  NAME 
!!   correct_electron_number
!!  USAGE
!! 
!!  PURPOSE
!!   Corrects the L matrix so that the electron number is correct
!!  INPUTS
!! 
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
!!  SOURCE
!!
  subroutine correct_electron_number( output_level, number_of_bands,inode, ionode)

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

    implicit none

    ! Passed variables
    integer :: inode, ionode, output_level

    real(double) ::  number_of_bands

    ! Local variables
    real(double) :: electrons_0, electrons, energy,  &
         step, step1, electrons2, g0, g1, recA, B, C, D, truestep, dne

    integer :: matTL, matphi2,matSphi,matSphi2
    logical :: done
    integer :: iter

    done = .false.
    iter = 0
    electrons_0 = two*number_of_bands
    matTL = allocate_temp_matrix(TLrange,0)
    matphi2 = allocate_temp_matrix(Lrange,0)
    matSphi = allocate_temp_matrix(Lrange,0)
    matSphi2 = allocate_temp_matrix(Lrange,0)
    ! get electron number and gradient
    call matrix_transpose(matT, matTtran)
   do while(.NOT.done.AND.(iter<20)) ! Was 20 !
    call LNV_matrix_multiply( electrons,energy, dontK, dontM1, dontM2, dontM3, dontM4, dophi, dontE,0,0,0,matphi)
    call matrix_product(matT,matphi,matTL,mult(T_L_TL))
    call matrix_product(matTL, matTtran,matSphi,mult(TL_T_L))
    if (inode .eq. ionode.and.output_level>=2) write(io_lun,1) electrons

    g0 = matrix_product_trace(matSphi, matphi)
    ! initial guess is linear correction...
    step1 = ( electrons_0 - electrons) / g0
    step = 0.1_double
    !step = step1
    if (inode.eq.ionode.and.output_level>=2) write(io_lun,*) 'g0, step1 are ',g0,step1
    ! if we are within 0.1% of the correct number, linear will do.
    !if (abs((electrons_0 - electrons)/electrons_0)<1e-6_double) then
    if (abs(electrons_0 - electrons)<1e-9_double) then
       call matrix_sum(one,matL,step1,matphi)
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
       call matrix_sum(one,matL,step1,matSphi)
       ! now get new electron number
       call LNV_matrix_multiply(electrons2,energy, &
            dontK, dontM1, dontM2, dontM3, dontM4, dophi, dontE,0,0,0,matphi2)
       call matrix_product(matT,matphi2,matTL,mult(T_L_TL))
       call matrix_product(matTL, matTtran,matSphi2,mult(TL_T_L))

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
  end subroutine correct_electron_number
!!***

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
