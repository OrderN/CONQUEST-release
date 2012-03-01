! $Id$
! -----------------------------------------------------------
! Module McWeeny
! -----------------------------------------------------------
! Code area 4: density matrix
! -----------------------------------------------------------

!!****h* Conquest/McWeeny
!!  NAME
!!   McWeeny
!!  PURPOSE
!!   To handle the McWeeny/Manolopoulos initialisation of density matrix
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   21/10/98
!!  MODIFICATION HISTORY
!!   11/05/01 DRB
!!    ROBODoc headers and better output
!!   08/06/2001 dave
!!    Added RCS Id and Log tags; used GenComms throughout
!!    Changed my_barrier to come from GenComms
!!   11:09, 04/02/2003 drb 
!!    Major change: stability fix in InitMcW to ensure mu lies between Hamiltonian limits
!!   10:09, 13/02/2006 drb 
!!    Removed all explicit references to data_ variables and rewrote in terms of new 
!!    matrix routines
!!   2006/11/21 16:53 dave
!!    Changed iprint statements
!!   2008/02/01 17:48 dave
!!    Changes for output to file not stdout
!!   2008/07/31 ast
!!    Added timers
!!   2012/03/01 L.Tong
!!    Added interface for McW_matrix_multiply
!!***
module McWeeny

  use global_module, ONLY: io_lun
  use timer_stdclocks_module, ONLY: start_timer, stop_timer, tmr_std_matrices

  implicit none

  !!****f* McWeeny/McW_matrix_multiply
  !! PURPOSE
  !!   Interface for McW_matrix_multiply_nondsw and
  !!   McW_matrix_multiply_ds
  !! AUTHOR
  !!   L.Tong
  !! CREATION DATE
  !!   2012/03/01
  !! MODIFICATION HISTORY
  !! SOURCE
  !!
  interface McW_matrix_multiply
     module procedure McW_matrix_multiply_nonds
     module procedure McW_matrix_multiply_ds
  end interface McW_matrix_multiply
  !!*****

contains

!!****f* McWeeny/McWMin *
!!
!!  NAME
!!   McWMin
!!  USAGE
!!
!!  PURPOSE
!!   Performs iterative minimisation using McWeeny algorithm
!!  INPUTS
!!
!!
!!  USES
!!   datatypes, numbers, logicals, maxima_module, matrix_data
!!   mult_module, multiply_module, GenBlas
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   21/10/98
!!  MODIFICATION HISTORY
!!   11/05/01 DRB
!!    ROBODoc headers and better output also nsf from maxima_module
!!   08/06/2001 dave
!!    Tidied and indented
!!   10:38, 16/05/2007 drb
!!    Moved matrix_sum to stop error at end
!!   2011/07/25  L.Tong
!!     Implemented spin polarisation
!!     added dependence on matL_dn, matH_dn, matLS, matSL, matL_dnS
!!     and matSL_dn
!!  SOURCE
!!
  subroutine McWMin (n_L_iterations, deltaE)

    use datatypes
    use numbers
    use matrix_data, ONLY: Lrange
    use mult_module, ONLY: allocate_temp_matrix, matrix_sum, &
         free_temp_matrix, matL, matH, matS, matT, symmetrise_L, &
         matL_dn, matH_dn, matLS, matSL, matL_dnS, matSL_dn
    use global_module, ONLY: iprint_DM, IPRINT_TIME_THRES1, &
         IPRINT_TIME_THRES2, flag_spin_polarisation, &
         flag_fix_spin_population
    use GenComms, ONLY: inode, ionode
    use timer_module, ONLY: cq_timer, start_timer, stop_print_timer, &
         WITH_LEVEL

    implicit none

    ! Passed variables
    integer :: n_L_iterations
    real(double) :: deltaE

    ! Local variables
    integer :: n_iterations
    integer :: matRhoNew, mat_oldL
    integer :: matRhoNew_dn, mat_oldL_dn
    real(double) :: cn, oldE, omega1, c_old
    real(double) :: cn_up, cn_dn, oldE_up, oldE_dn, omega1_up, &
         omega1_dn, c_old_up, c_old_dn
    type(cq_timer) :: tmr_l_tmp1,tmr_l_tmp2

    call start_timer(tmr_l_tmp1,WITH_LEVEL)
    oldE = 1.0e30_double
    if (flag_spin_polarisation .and. flag_fix_spin_population) then
       oldE_up = 1.0e30_double
       oldE_dn = 1.0e30_double
    end if
    mat_oldL = allocate_temp_matrix(Lrange,0)
    matRhoNew = allocate_temp_matrix(Lrange,0)
    call matrix_sum (zero, mat_oldL, one, matL)
    if (flag_spin_polarisation) then
       mat_oldL_dn = allocate_temp_matrix (Lrange,0)
       matRhoNew_dn = allocate_temp_matrix (Lrange, 0)
       call matrix_sum (zero, mat_oldL_dn, one, matL_dn)
       if (flag_fix_spin_population) then ! do each spin channel separately
          do n_iterations=1,n_L_iterations
             call start_timer(tmr_l_tmp2,WITH_LEVEL)
             if (inode==ionode.and.iprint_DM>=1) write(io_lun,2) n_iterations
             ! find new rho and energy
             call McW_matrix_multiply (matH, matS, matL, matT, matLS, &
                  matSL, matRhoNew, omega1_up, cn_up)
             call McW_matrix_multiply (matH_dn, matS, matL_dn, matT, &
                  matL_dnS, matSL_dn, matRhoNew_dn, omega1_dn, cn_dn)
             ! If new energy is lower, copy new rho into old rho and repeat
             if ((omega1_up > oldE_up .or. omega1_dn > oldE_dn) .or. &
                  ((abs (omega1_up - oldE_up) < 1.0e-08) .or. &
                  (abs (omega1_dn - oldE_dn) < 1.0e-08)) .or. &
                  (cn_up > one .or. cn_dn > one) .or. &
                  (cn_up < zero .or. cn_dn < zero)) then
                if(inode==ionode.AND.iprint_DM>0) &
                     write(io_lun, 6) c_old_up, c_old_dn, oldE_up, &
                     oldE_dn, cn_up, cn_dn, omega1_up, omega1_dn
                call stop_print_timer(tmr_l_tmp2, "a McWeeny iteration",&
                     IPRINT_TIME_THRES2)
                exit
             end if
             deltaE = (omega1_up - oldE_up) + (omega1_dn - oldE_dn)
             call matrix_sum (zero, matL, one, matRhoNew)
             call matrix_sum (zero, matL_dn, one, matRhoNew_dn)
             oldE_up = omega1_up
             oldE_dn = omega1_dn
             c_old_up = cn_up
             c_old_dn = cn_dn
             if(inode==ionode.AND.iprint_DM>=0) &
                  write(io_lun, 7) oldE_up, oldE_dn, c_old_up, c_old_dn
             call symmetrise_L (spin=1)
             call symmetrise_L (spin=2)
             call stop_print_timer (tmr_l_tmp2, &
                  "a McWeeny iteration", IPRINT_TIME_THRES2)
          end do
       else ! do each spin channel as a direct sum
          do n_iterations=1,n_L_iterations
             call start_timer(tmr_l_tmp2,WITH_LEVEL)
             if (inode==ionode.and.iprint_DM>=1) write(io_lun,2) n_iterations
             ! find new rho and energy
             call McW_matrix_multiply (matH, matH_dn, matS, &
                  matL, matL_dn, matT, matRhoNew, matRhoNew_dn, &
                  omega1, cn)
             ! If new energy is lower, copy new rho into old rho and repeat
             ! If the spin populations are allowed to vary there will
             ! be only one omega1 and cn etc for the entire direct sum
            if (omega1 > oldE .or. &
                 (abs (omega1 - oldE) < 1.0e-8) .or. &
                 cn > one .or. cn < zero) then
               if(inode==ionode.AND.iprint_DM>0) &
                    write(io_lun,4) c_old, oldE, cn, omega1
               call stop_print_timer(tmr_l_tmp2, "a McWeeny iteration",&
                    IPRINT_TIME_THRES2)
               exit
             end if
             deltaE = omega1 - oldE
             call matrix_sum (zero, matL, one, matRhoNew)
             call matrix_sum (zero, matL_dn, one, matRhoNew_dn)
             oldE = omega1
             c_old = cn
             if(inode==ionode.AND.iprint_DM>=0) write(io_lun,5) oldE, c_old
             call symmetrise_L (spin=1)
             call symmetrise_L (spin=2)
             call stop_print_timer (tmr_l_tmp2, &
                  "a McWeeny iteration", IPRINT_TIME_THRES2)
          end do
       end if ! if (flag_fix_spin_population)
    else  ! spin non-polarised calculation
       do n_iterations= 1, n_L_iterations
          call start_timer (tmr_l_tmp2, WITH_LEVEL)
          if (inode==ionode .and. iprint_DM>=1) write (io_lun, 2) n_iterations
          ! find new rho and energy
          call McW_matrix_multiply (matH, matS, matL, matT, matLS, &
               matSL, matRhoNew, omega1, cn)
          ! If new energy is lower, copy new rho into old rho and repeat
          if(omega1.GT.oldE.OR.&
               (abs(omega1-oldE)<1.0e-8).OR.&
               cn.GT.one.OR.cn.LT.zero) then
             if(inode==ionode.AND.iprint_DM>0) &
                  write(io_lun,4) c_old, oldE, cn, omega1
             call stop_print_timer(tmr_l_tmp2, "a McWeeny iteration", &
                  IPRINT_TIME_THRES2)
             exit
          endif
          deltaE = omega1 - oldE
          call matrix_sum (zero, matL, one, matRhoNew)
          oldE = omega1
          c_old = cn
          if(inode==ionode.AND.iprint_DM>=0) write(io_lun,5) oldE, c_old
          call symmetrise_L()
          call stop_print_timer (tmr_l_tmp2,&
               "a McWeeny iteration",IPRINT_TIME_THRES2)
       end do
    end if
    if(n_iterations<3) call matrix_sum(zero,matL,one,mat_oldL)
    if (flag_spin_polarisation) then
       call free_temp_matrix (matRhoNew_dn)
       call free_temp_matrix (mat_oldL_dn)
    end if
    call free_temp_matrix(matRhoNew)
    call free_temp_matrix(mat_oldL)
    if (inode.eq.ionode.and.iprint_DM.ge.1) write(io_lun,3) n_iterations, omega1
    call stop_print_timer(tmr_l_tmp1,"MCWEENY",IPRINT_TIME_THRES1)
2   format(/,20x,'McWeeny L iteration:',i5)
3   format(/,20x,'Functional value reached after ',i5,' L iterations: ', &
         /,20x,' Omega: ', f15.7)
4   format(2x,'Rounding error found: ', 4f15.6)
5   format(3x,'Energy: ',f15.6,' Mid-point: ',f15.6)
6   format(2x,'Rounding error found: ', 8f15.6)
7   format(3x,'Energy (up and down): ', 2f15.6,' Mid-point (up and down): ', 2f15.6)
    return
  end subroutine McWMin
!!***

!!****f* McWeeny/InitMcW *
!!
!!  NAME 
!!   InitMcW 
!!  USAGE
!! 
!!  PURPOSE
!!   Initialises L to be a function of H with eigenvalues between 0 and 1
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   End 10/98
!!  MODIFICATION HISTORY
!!   11/05/01 DRB
!!    ROBODoc headers and better output also nsf from maxima_module
!!   08/06/2001 dave
!!    Tidying and added GenComms for gsum
!!    Added my_barrier to use GenComms statement
!!   15:25, 10/12/2002 drb 
!!    Bug fix - changed location of output for mubar and added gmin and gmax calls on hmin and hmax
!!   11:08, 04/02/2003 drb 
!!    Implemented quick fix for bracketing mu with Hamiltonian limits (required for stability)
!!   12:30, 2004/06/09 dave
!!    Fixed assumption that T and S have same range
!!   2004/10/29 drb
!!    Removed gsum on SX: MAJOR bug found by TM which was messing up initialisation;
!!    also removed loops to change hmin and hmax (unnecessary, I think)
!!   2011/07/25 L.Tong 
!!    Adding Spin polarisation
!!   Friday, 2011/08/05 L.Tong 
!!    Removed the obsolete number_of_bands parameter, and the total
!!    electron numbers are obtained directly from the global_module
!!   2012/01/30 L.Tong
!!    Rewrite for the variable spin case, works for the whole direct sum.
!!  SOURCE
!!
  subroutine InitMcW (inode, ionode)

    use datatypes
    use numbers
    use global_module, ONLY: iprint_DM, ni_in_cell, &
         IPRINT_TIME_THRES1, flag_spin_polarisation, &
         flag_fix_spin_population, ne_in_cell, ne_up_in_cell, &
         ne_dn_in_cell
    use species_module, ONLY: species, nsf_species
    use mult_module, ONLY: allocate_temp_matrix, matH, matS, matT, &
         matL, mult, matrix_trace, matrix_product_trace, &
         matrix_product, T_S_TS, free_temp_matrix, matH_dn, matL_dn
    use matrix_data, ONLY: Lrange, Srange, TSrange
    use GenBlas
    use GenComms, ONLY: gsum, my_barrier, gmin, gmax
    use timer_module, ONLY: cq_timer,start_timer, stop_print_timer, &
         WITH_LEVEL

    implicit none

    ! Passed variables
    ! real(double) :: number_of_bands 
    integer :: inode, ionode

    ! Local variables
    real(double) :: n_e, n_o, hmax, hmin, SX, SXHX, SXH_upX, SXH_dnX
    real(double) :: A, mu1, mu2, mubar, lambda
    real(double) :: ne_up, mu1_up, mu2_up, mubar_up, lambda_up, hmin_up, hmax_up
    real(double) :: ne_dn, mu1_dn, mu2_dn, mubar_dn, lambda_dn, hmin_dn, hmax_dn
    integer :: matXHX, mat_temp, matTS
    integer :: matXH_dnX, mat_temp_dn
    integer :: length, i
    type(cq_timer) :: tmr_l_tmp1

    call start_timer(tmr_l_tmp1,WITH_LEVEL)
    if(inode==ionode.AND.iprint_DM>1) write(io_lun,1)
    ! We must first initialise rho
    ! n_e is the number of electrons per spin channel
    if (flag_spin_polarisation) then
       if (flag_fix_spin_population) then
          ne_up = ne_up_in_cell
          ne_dn = ne_dn_in_cell
       else ! variable spin
          n_e = ne_in_cell
       end if
    else ! no spin
       n_e = half * ne_in_cell
    end if

    ! n_o is number of support functions 
    n_o = zero
    do i=1, ni_in_cell
       n_o = n_o + nsf_species(species(i))
    end do
    
    ! we are working with the whole direct sum for variable spin case
    ! so we need to have twice the size of basis set
    if (flag_spin_polarisation .and. &
         (.not. flag_fix_spin_population)) &
         n_o = two * n_o

    if (flag_fix_spin_population) then
       if (inode == ionode .AND. iprint_DM > 1) &
            write (io_lun, 8) ne_up, ne_dn, n_o
    else
       if (inode == ionode .AND. iprint_DM > 1) &
            write (io_lun, 2) n_e, n_o
    end if
    
    matXHX = allocate_temp_matrix(Lrange,0)
    mat_temp = allocate_temp_matrix(Srange,0)
    matTS = allocate_temp_matrix(TSrange,0)
    if (flag_spin_polarisation) then
       matXH_dnX = allocate_temp_matrix(Lrange,0)
       mat_temp_dn = allocate_temp_matrix(Srange,0)
    end if
    call McWXHX (matH, matT, matXHX, mat_temp)
    call matrix_product (matT, matS, matTS, mult(T_S_TS))

    SX = matrix_trace (matTS)

    ! for variable spin, we are working with direct sum, so need to
    ! take the sum of spin up and dn components
    if (flag_spin_polarisation .and. &
         (.not. flag_fix_spin_population)) SX = two * SX

    if (inode==ionode.AND.iprint_DM>1) write (io_lun,*) 'SX is ', SX
    if (inode==ionode.AND.iprint_DM>1) write (io_lun,*) 'SX/n_o is ', SX/n_o

    if (flag_spin_polarisation) then
       call McWXHX (matH_dn, matT, matXH_dnX, mat_temp_dn)
       SXH_upX = matrix_product_trace (matS, mat_temp)
       SXH_dnX = matrix_product_trace (matS, mat_temp_dn)
       if (inode==ionode.AND.iprint_DM>1) &
            write (io_lun,*) 'SXH_upX, SXH_dnX are ', SXH_upX, SXH_dnX
       if (.not. flag_fix_spin_population) then
          ! total trace = sum of traces
          SXHX = SXH_upX + SXH_dnX
          if (inode==ionode.AND.iprint_DM>1) write (io_lun,*) 'SXHX is ', SXHX
       end if
    else
       SXHX = matrix_product_trace (matS,mat_temp)
       if (inode==ionode.AND.iprint_DM>1) write (io_lun,*) 'SXHX is ', SXHX
    end if
    call my_barrier()
    ! Get Gershgorin limits on H
    if (flag_spin_polarisation) then
       call GetHLimits (hmin_up, hmax_up, spin=1)
       call GetHLimits (hmin_dn, hmax_dn, spin=2)
       call gmin (hmin_up)
       call gmin (hmin_dn)
       call gmax (hmax_up)
       call gmax (hmax_dn)
       if (inode==ionode.AND.iprint_DM>1) then
          write (io_lun, 4) hmin_up, hmax_up
          write (io_lun, 5) hmin_dn, hmax_dn
       end if
       if (.not. flag_fix_spin_population) then
          ! we need to find the min and max across the spin channels
          ! since we are treating the overall matrix as direct sum of the
          ! spin parts and store them in hmin and hmax
          hmin = min (hmin_up, hmin_dn)
          hmax = max (hmax_up, hmax_dn)
          if (inode==ionode.AND.iprint_DM>1) write (io_lun, 6) hmin, hmax
       end if
    else
       call GetHLimits (hmin, hmax)
       call gmin (hmin)
       call gmax (hmax)
       if (inode==ionode.AND.iprint_DM>1) write (io_lun, 6) hmin, hmax
    end if

    A = (one - SX/n_o)
    if (flag_fix_spin_population) then
       mu1_up = SXH_upX / n_o + hmax_up * A
       mu1_dn = SXH_dnX / n_o + hmax_dn * A
       mu2_up = (hmin_up * A * ne_up / (n_o - ne_up) - SXH_upX / n_o) &
            / (n_o * A / (n_o - ne_up) - one)
       mu2_dn = (hmin_dn * A * ne_dn / (n_o - ne_dn) - SXH_dnX / n_o) & 
            / (n_o * A / (n_o - ne_dn) - one)
       if(inode==ionode) write(io_lun,*) 'Mu1_up, Mu_dn, Mu2_up, Mu2_dn: ', &
            mu1_up, mu1_dn, mu2_up, mu2_dn
       ! spin up channel
       if ((mu1_up < hmax_up .AND. mu1_up > hmin_up) &
            .AND. (abs (ne_up / (hmax_up - mu1_up)) .LT. &
            abs ((n_o - ne_up) / (mu2_up - hmin_up)))) then
          mubar_up = mu1_up
          lambda_up = ne_up / (hmax_up - mu1_up)
       else
          mubar_up = mu2_up
          lambda_up = (n_o - ne_up) / (mu2_up - hmin_up)
       endif
       ! spin down channel
       if ((mu1_dn < hmax_dn .AND. mu1_dn > hmin_dn) &
            .AND. (abs (ne_dn / (hmax_dn - mu1_dn)) .LT. &
            abs ((n_o - ne_dn) / (mu2_dn - hmin_dn)))) then
          mubar_dn = mu1_dn
          lambda_dn = ne_dn / (hmax_dn - mu1_dn)
       else
          mubar_dn = mu2_dn
          lambda_dn = (n_o - ne_dn) / (mu2_dn - hmin_dn)
       endif
       if(inode==ionode.AND.iprint_DM>1) write(io_lun,*) &
            'mubar_up, mubar_dn are ', mubar_up, mubar_dn
       if(inode==ionode.AND.iprint_DM>1) write(io_lun,*) &
            'lambda_up and lambda_dn are ', lambda_up, lambda_dn
    else ! if (flag_fix_spin_population) then
       mu1 = SXHX / n_o + hmax * A
       mu2 = (hmin * A * n_e / (n_o - n_e) - SXHX / n_o) / &
            (n_o * A / (n_o - n_e) - one)
       if (inode == ionode) write(io_lun,*) 'Mu1,2: ', mu1, mu2
       ! DRB 2004/09/16 Fixing typo ?
       if ((mu1 < hmax .AND. mu1 > hmin) .AND. &
            (abs(n_e / (hmax - mu1)) .LT. &
            abs((n_o - n_e) / (mu2 - hmin)))) then
          mubar = mu1
          lambda = n_e / (hmax - mu1)
       else
          mubar = mu2
          lambda = (n_o - n_e) / (mu2 - hmin)
       endif
       if (inode == ionode .AND. iprint_DM > 1) &
            write (io_lun, 7) mubar
       if (inode == ionode .AND. iprint_DM > 1) &
            write (io_lun, *) 'lambda is ', lambda
    end if ! if (flag_fix_spin_population) then
    
    ! Calculate L0. 
    if (flag_spin_polarisation) then
       if (flag_fix_spin_population) then
          call McWRho0 (matL, matT, matXHX, mubar_up, lambda_up, ne_up, n_o)
          call McWRho0 (matL_dn, matT, matXH_dnX, mubar_dn, lambda_dn, ne_dn, n_o)
       else
          call McWRho0 (matL, matT, matXHX, mubar, lambda, n_e, n_o)
          call McWRho0 (matL_dn, matT, matXH_dnX, mubar, lambda, n_e, n_o)
       end if
    else
       call McWRho0 (matL, matT, matXHX, mubar, lambda, n_e, n_o)
    end if

    ! Free the temporary matrices
    if (flag_spin_polarisation) then
       call free_temp_matrix(mat_temp_dn)
       call free_temp_matrix(matXH_dnX)
    end if
    call free_temp_matrix(matTS)
    call free_temp_matrix(mat_temp)
    call free_temp_matrix(matXHX)
    call stop_print_timer(tmr_l_tmp1, "McWeeny initialisation", IPRINT_TIME_THRES1)

1   format(1x,'Welcome to InitMcW')
2   format(2x,'Electrons: ',f15.6,' Orbitals: ',f15.6)
8   format(2x,'Electrons_up: ',f15.6,' Electrons_dn: ',f15.6,' Orbitals: ',f15.6)
3   format(2x,'Electrons_up: ',f15.6,' Electrons_dn: ',f15.6)
4   format(2x,'Minimum and maximum limits on H_up are ',2f15.6)
5   format(2x,'Minimum and maximum limits on H_dn are ',2f15.6)
6   format(2x,'Minimum and maximum limits on H are ',2f15.6)
7   format(2x,'Mubar is ',f15.6)

    return
  end subroutine InitMcW
!!***

!!****f* McWeeny/McW_matrix_multiply_nonds *
!!
!!  NAME
!!   McW_matrix_multiply_nonds
!!  USAGE
!!
!!  PURPOSE
!!   Performs one iteration of McWeeny iterative scheme
!!   Ref: PRB 58, 12704 (1998)
!!  INPUTS
!!
!!  USES
!!
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   21/10/98
!!  MODIFICATION HISTORY
!!   11/05/01 DRB
!!    ROBODoc headers and better output also nsf from maxima_module
!!   08/06/2001 dave
!!    Added GenComms for gsum
!!   2011/03/28 L.Tong
!!    Removed local variable declaration matL2, redundant
!!   2011/07/26 L.Tong
!!    Implemented spin polarisation. For this subroutine, the spins
!!    will be treated separately all together, and this subroutine
!!    should be called twice involking different spin The matH, matL
!!    and matRhoNew, and energy, c should be used accordingly for
!!    different spin channels.
!!    Decided for simplification, moving matLS and matSL refereces to
!!    be part of the subroutine parameters instead of module
!!    references. This tidies the code somewhat, in this case the
!!    usual optional spin parameter is not required.
!!    But still needs the flag_spin_polarisation flag to toggle the
!!    factor of two for total energy calculation
!!   2012/03/01 L.Tong
!!    Changed name to McW_matrix_multiply_nonds, and to be called from
!!    interface McW_matrix_multiply
!!  SOURCE
!!
  subroutine McW_matrix_multiply_nonds(matH, matS, matL, matInvS, &
       matLS, matSL, matRhoNew, energy, c)

    use datatypes
    use numbers
    use primary_module, ONLY : bundle
    use matrix_data, ONLY: LSLrange, Lrange, Hrange, Srange
    use mult_module, ONLY: matrix_product, matrix_sum, &
         matrix_product_trace, matrix_scale, matrix_transpose, &
         allocate_temp_matrix, free_temp_matrix, L_S_LS, LS_L_LSL, &
         LSL_SL_L, mult
    use global_module, ONLY: iprint_DM, flag_spin_polarisation
    use GenComms, ONLY: gsum, cq_abort, inode, ionode, my_barrier

    implicit none

    ! Passed variables
    integer :: matH, matS, matL, matInvS, matRhoNew
    integer :: matLS, matSL
    real( double ) energy, c

    ! Local Variables
    real( double ) c1,c2,cn
    integer :: matLSL, matLSLSL, mat_top, mat_bottom

    call matrix_product(matL, matS, matLS, mult(L_S_LS))
    ! Generate LSL, LSLSL
    matLSL = allocate_temp_matrix(LSLrange,0)
    matLSLSL = allocate_temp_matrix(Lrange,0)
    mat_top = allocate_temp_matrix(Srange,0)
    mat_bottom = allocate_temp_matrix(Srange,0)
    call my_barrier
    call matrix_product(matLS,matL,matLSL,mult( LS_L_LSL))
    call matrix_product(matLSL,matLS,matLSLSL,mult( LSL_SL_L))
    !----------------------------------------------------------------
    ! Shorten L, LSL and LSLSL to range S and take traces (can do it
    ! with sdot because S is symmetric)
    !----------------------------------------------------------------
    call matrix_sum(zero,mat_top,one,matLSL)
    call matrix_sum(one,mat_top,-one,matLSLSL)
    c1 = matrix_product_trace(matS,mat_top)
    if(inode==ionode.AND.iprint_DM>=2) write(io_lun,*) 'S.top is ',c1
    call matrix_sum(zero,mat_bottom,one,matL)
    c1 = matrix_product_trace(matS,mat_bottom)
    if(inode==ionode.AND.iprint_DM>=2) write(io_lun,*) 'N_e is ',c1
    call matrix_sum(one,mat_bottom,-one,matLSL)
    c1 = matrix_product_trace(matS,mat_top)
    c2 = matrix_product_trace(matS,mat_bottom)
    if(c2/=zero) then
       cn = c1/c2
    else
       call cq_abort('McW_matrix_multiply: c2 is zero')
    endif
    c = cn
    if(inode.eq.ionode.AND.iprint_DM>=2) write(io_lun,*) 'c, c1,c2 are ',cn,c1,c2
    ! Shorten LSL and LSLSL to range L
    call matrix_scale(zero,matRhoNew)
    call matrix_sum(zero,matRhoNew,-one,matLSLSL)
    call matrix_sum(one,matRhoNew,one+cn,matLSL)
    if(cn<0.5_double) then
       call matrix_sum(one,matRhoNew,(one - two*cn),matL)
       call matrix_scale((one/(one - cn)),matRhoNew)
    else
       call matrix_scale(one/cn,matRhoNew)
    endif
    ! Idempotency checks
    !%%! data_L3 = zero
    !%%! call axpy(NSF*NSF*mat(1,Lrange)%length,one,dat_L,1,data_L3,1)
    !%%! ! LSLSL -> 2LSLSL - 3LSL
    !%%! call matrix_add(2.0_double, data_LSLSL, mat(:,Lrange), &
    !%%!      -3.0_double, data_LSL, mat(:,LSLrange), bundle, inode)
    !%%! ! Lnew -> Lnew+(2LSLSL-3LSL) [= Lnew - McW(Lold)]
    !%%! call axpy(NSF*NSF*mat(1,Lrange)%length,one,data_LSLSL,1,data_L3,1)
    !%%! call mult_wrap(inode-1,data_L3,dat_S,data_LS,mult(L_S_LS))
    !%%! call mult_wrap(inode-1,dat_S,data_SL,data_XHX,mult(S_LS_L))
    !%%! cn = dot(NSF*NSF*mat(1,Lrange)%length,data_SLS2, 1, data_L3, 1)
    !%%! call gsum(cn)
    !%%! if(inode.eq.ionode.AND.iprint_DM>=2) write(io_lun,*) 'Idempot monitor is ',cn/256.0_double
    ! Shorten rhonew to range H into temp
    call free_temp_matrix(mat_bottom)
    mat_bottom = allocate_temp_matrix(Hrange,0)
    call matrix_sum(zero,mat_bottom,one,matRhoNew)
    if (flag_spin_polarisation) then
       energy = matrix_product_trace(mat_bottom,matH)
    else
       energy = two*matrix_product_trace(mat_bottom,matH)
    end if
    if(inode==ionode.AND.iprint_DM>=3) write(io_lun,*) 'energy is ',energy
    call matrix_sum(zero,mat_top,one,matRhoNew)
    c1 = matrix_product_trace(mat_top,matS)
    if(inode==ionode.AND.iprint_DM>=3) write(io_lun,*) 'N_e(2) is ',c1
    ! note that N_e is printed without factor of two here even for spin un-polarised results.
    call free_temp_matrix(mat_bottom)
    call free_temp_matrix(mat_top)
    call free_temp_matrix(matLSLSL)
    call free_temp_matrix(matLSL)
    if(abs(c2)<1e-6) c = -0.01_double
  End Subroutine McW_matrix_multiply_nonds
!!***


!!****f* McWeeny/McW_matrix_multiply_ds *
!! PURPOSE
!!   Performs one iteration of McWeeny iterative scheme treating the
!!   overall matrices L, H as direct sums of their spin
!!   components. This allows population migration from one spin
!!   channel to another
!!   Ref: PRB 58, 12704 (1998)
!! INPUTS
!! OUTPUT
!! RETURN VALUE
!! AUTHOR
!!   L.Tong
!! CREATION DATE
!!   2011/07/26
!! MODIFICATION HISTORY
!!   Tuesday, 2011/08/02 L.Tong
!!     Changed the order of input parameters, replaced energy_up and
!!     energy_down with energy_total, and c with c_total
!!   2012/03/01 L.Tong
!!     Renamed to McW_matrix_multiply_ds, to be called by interface
!!     McW_matrix_multiply
!! SOURCE
!!
  subroutine McW_matrix_multiply_ds (matH, matH_down, matS, matL, matL_down, &
       matInvS, matRhoNew, matRhoNew_down, energy_total, c_total)

    use datatypes
    use numbers
    use primary_module, ONLY : bundle
    use matrix_data, ONLY: LSLrange, Lrange, Hrange, Srange
    use mult_module, ONLY: matrix_product, matrix_sum, &
         matrix_product_trace, matrix_scale, matrix_transpose, &
         allocate_temp_matrix, free_temp_matrix, matLS, matSL, L_S_LS,&
         LS_L_LSL, LSL_SL_L, mult, matL_dnS, matSL_dn
    use global_module, ONLY: iprint_DM, flag_spin_polarisation
    use GenComms, ONLY: gsum, cq_abort, inode, ionode, my_barrier

    implicit none

    ! Passed variables
    integer :: matH, matS, matL, matInvS, matRhoNew
    integer :: matH_down, matL_down, matRhoNew_down
    real(double) :: energy_total, c_total

    ! Local Variables
    real(double) :: c1,c2,c3,cn
    integer :: matLSL, matLSLSL, mat_top, mat_bottom
    integer :: matL_dnSL_dn, matL_dnSL_dnSL_dn

    call matrix_product (matL, matS, matLS, mult(L_S_LS))
    call matrix_product (matL_down, matS, matL_dnS, mult(L_S_LS))
    ! Generate LSL, LSLSL
    matLSL = allocate_temp_matrix (LSLrange, 0)
    matLSLSL = allocate_temp_matrix (Lrange, 0)
    matL_dnSL_dn = allocate_temp_matrix (LSLrange, 0)
    matL_dnSL_dnSL_dn = allocate_temp_matrix (Lrange,0)
    mat_top = allocate_temp_matrix (Srange, 0)
    mat_bottom = allocate_temp_matrix (Srange, 0)
    call my_barrier
    call matrix_product (matLS, matL, matLSL, mult(LS_L_LSL))
    call matrix_product (matLSL, matLS, matLSLSL, mult( LSL_SL_L))
    call matrix_product (matL_dnS, matL_down, matL_dnSL_dn, mult( LS_L_LSL))
    call matrix_product (matL_dnSL_dn, matL_dnS, matL_dnSL_dnSL_dn, mult( LSL_SL_L))
    !----------------------------------------------------------------
    ! Shorten L, LSL and LSLSL to range S and take traces (can do it
    ! with sdot because S is symmetric)
    !----------------------------------------------------------------
    ! spin up contribution
    call matrix_sum (zero, mat_top, one, matLSL)
    call matrix_sum (one, mat_top, -one, matLSLSL)
    c1 = matrix_product_trace (matS, mat_top)
    if (inode==ionode .AND. iprint_DM>=2) write (io_lun,*) 'S.top_up is ', c1
    call matrix_sum (zero, mat_bottom, one, matL)
    c1 = matrix_product_trace (matS, mat_bottom)
    if (inode==ionode .AND. iprint_DM>=2) write(io_lun,*) 'N_e_up is ', c1
    ! get the real mat_bottom for spin up
    call matrix_sum (one, mat_bottom, -one, matLSL)
    ! real c1 and c2 calculated
    c1 = matrix_product_trace (matS, mat_top)
    c2 = matrix_product_trace (matS, mat_bottom)
    ! spin down contribution
    call matrix_sum (zero, mat_top, one, matL_dnSL_dn)
    call matrix_sum (one, mat_top, -one, matL_dnSL_dnSL_dn)
    ! c1 and c2 cannot be used as temporary here, so use a new local variable
    c3 = matrix_product_trace (matS, mat_top)
    if (inode==ionode.AND.iprint_DM>=2) write (io_lun,*) 'S.top_dn is ',c3
    call matrix_sum (zero, mat_bottom, one, matL_down)
    c3 = matrix_product_trace (matS, mat_bottom)
    if (inode==ionode.AND.iprint_DM>=2) write (io_lun,*) 'N_e_dn is ',c3
    ! get the real mat_bottom for spin down
    call matrix_sum (one, mat_bottom, -one, matL_dnSL_dn)
    ! c1 and c2 accumulated
    c1 = c1 + matrix_product_trace (matS, mat_top)
    c2 = c2 + matrix_product_trace (matS, mat_bottom)
    if (c2/=zero) then
       cn = c1/c2
    else
       call cq_abort ('McW_matrix_multiply: c2 is zero')
    endif
    c_total = cn
    if (inode.eq.ionode.AND.iprint_DM>=2) &
         write (io_lun,*) 'c_total, c1,c2 are ', cn, c1, c2
    ! c1 and c2 is now free to be used as temporary variables again
    ! Shorten LSL and LSLSL to range L
    call matrix_scale (zero, matRhoNew)
    call matrix_sum (zero, matRhoNew, -one, matLSLSL)
    call matrix_sum (one, matRhoNew, one+cn, matLSL)
    if (cn<0.5_double) then
       call matrix_sum (one, matRhoNew,(one - two*cn),matL)
       call matrix_scale ((one/(one - cn)),matRhoNew)
    else
       call matrix_scale (one/cn, matRhoNew)
    endif
    ! do the second spin channel
    call matrix_scale (zero, matRhoNew_down)
    call matrix_sum (zero, matRhoNew_down, -one, matL_dnSL_dnSL_dn)
    call matrix_sum (one, matRhoNew_down, one+cn, matL_dnSL_dn)
    if (cn<0.5_double) then
       call matrix_sum (one, matRhoNew_down,(one - two*cn),matL_down)
       call matrix_scale ((one/(one - cn)),matRhoNew_down)
    else
       call matrix_scale (one/cn, matRhoNew_down)
    endif
    call free_temp_matrix (mat_bottom)
    mat_bottom = allocate_temp_matrix (Hrange,0)
    call matrix_sum (zero, mat_bottom, one, matRhoNew)
    energy_total = matrix_product_trace (mat_bottom, matH)
    call matrix_sum (zero, mat_bottom, one, matRhoNew_down)
    energy_total = energy_total + matrix_product_trace (mat_bottom, matH_down)
    if (inode==ionode.AND.iprint_DM>=3) write (io_lun,*) 'energy_total ', energy_total
    call matrix_sum (zero, mat_top, one, matRhoNew)
    c1 = matrix_product_trace (mat_top, matS)
    call matrix_sum (zero, mat_top, one, matRhoNew_down)
    c2 = matrix_product_trace (mat_top, matS)
    if (inode==ionode.AND.iprint_DM>=3) write (io_lun,*) 'N_e_up, N_e_dn are ', c1, c2
    ! free temp matrices
    call free_temp_matrix (mat_bottom)
    call free_temp_matrix (mat_top)
    call free_temp_matrix (matL_dnSL_dnSL_dn)
    call free_temp_matrix (matL_dnSL_dn)
    call free_temp_matrix (matLSLSL)
    call free_temp_matrix (matLSL)
  End Subroutine McW_matrix_multiply_ds
!!***

!!****f* MacWeeny/McWXHX
!! PURPOSE
!!  sbrt McWXHX: Pre- and post- multiplies H by S^-1, and returns the result
!!  to both L range (bab) and S range (tm)
!! INPUTS
!! OUTPUT
!! RETURN VALUE
!! AUTHOR
!!   drb
!! CREATION DATE
!!
!! MODIFICATION HISTORY
!! SOURCE
!!
  subroutine McWXHX(matA, matB, matBAB, mat_temp)

    use datatypes
    use numbers
    use matrix_data, ONLY: THrange
    use mult_module, ONLY: matrix_product, matrix_sum, T_H_TH, TH_T_L, &
         allocate_temp_matrix, mult, free_temp_matrix, matrix_trace
    use GenComms, ONLY : inode, ionode
    use global_module, ONLY: iprint_DM

    implicit none

    ! Passed variables
    integer :: matA, matB, matBAB, mat_temp

    ! Local variable
    integer :: matBA
    real(double) :: xx

    matBA = allocate_temp_matrix(THrange,0)
    call matrix_product(matB, matA, matBA, mult(T_H_TH))
    if (iprint_DM >= 2) xx = matrix_trace(matBA)
    if (inode .eq. ionode .AND. iprint_DM >= 2) &
         write (io_lun, *) 'Trace of BA: ', xx
    call matrix_product(matBA, matB, matBAB, mult(TH_T_L))
    if (iprint_DM >= 2) xx = matrix_trace(matBAB)
    if (inode .eq. ionode .AND. iprint_DM >= 2) &
         write (io_lun, *) 'Trace of BAB: ', xx
    call matrix_sum(zero, mat_temp, one, matBAB)
    if (iprint_DM >= 2) xx = matrix_trace(mat_temp)
    if (inode .eq. ionode .AND. iprint_DM >= 2) &
         write (io_lun, *) 'Trace of BAB: ', xx
    call free_temp_matrix(matBA)
    return
  end Subroutine McWXHX
!!*****

!!****f*
!! PURPOSE
!! INPUTS
!! OUTPUT
!! RETURN VALUE
!! AUTHOR
!!   drb
!! CREATION DATE
!! MODIFICATION HISTORY
!! SOURCE
!!
  subroutine McWRho0(matA, matB, matC, m, l, n_e, n_o)

    use datatypes
    use numbers
    use GenComms, ONLY : inode, ionode
    use mult_module, ONLY : matrix_sum
    use global_module, ONLY: iprint_DM

    implicit none

    integer :: matA, matB, matC
    Real( double ) :: m, l, n_e, n_o, tmp

    tmp = (l*m+n_e)/n_o
    if (inode .eq. ionode .AND. iprint_DM >= 2) &
         write (io_lun, *) 'l,m,n_e,n_o and tmp are ', &
         l, m, n_e, n_o, tmp
    call matrix_sum(zero, matA, tmp, matB)
    tmp = -l/n_o
    if (inode .eq. ionode .AND. iprint_DM >= 2) &
         write (io_lun, *) 'tmp is ', tmp
    call matrix_sum(one, matA, tmp, matC)
    return
  end Subroutine McWRho0
!!*****

!!****f*
!! PURPOSE
!!   Finds the Gershgorin limits on the eigenvalues of the Hamiltonian
!! INPUTS
!! OUTPUT
!! RETURN VALUE
!! AUTHOR
!!   drb
!! CREATION DATE
!! MODIFICATION HISTORY
!!   2011/07/25 L.Tong
!!     Implemented Spin Polarisation support: optional variable spin, that
!!     makes the subroutine to calculate spin up of spin down components.
!!   Wednesday, 2011/08/03 L.Tong
!!     Fixed the bug caused by the conditional if (present (spin)
!!     .and. spin == 2), on some compilers which does not short circuit
!!     and evaluates both sides of .and. and when the subroutine is
!!     called without the optional spin parameter, then spin == 2 gives
!!     a segmentation fault. Fixed this by introducing logical variable
!!     do_spin_down, we work out this variable correctly first, and
!!     then apply conditional on it in place of the previous present
!!     (spin) .and. spin == 2 conditional. Also moved the conditional
!!     outside of the inner loops to improve efficiency
!! SOURCE
!!
  subroutine GetHLimits(MinH, MaxH, spin)

    use datatypes
    use matrix_module, ONLY: matrix, matrix_halo
    use matrix_data, ONLY: mat, Hrange, halo
    use cover_module, ONLY: BCS_parts
    use primary_module, ONLY: bundle
    use mult_module, ONLY: matH, return_matrix_value, matrix_pos, matH_dn
    use numbers
    use GenComms,ONLY: cq_abort
    use maxima_module, ONLY: maxnsf

    implicit none

    ! Passed variables
    real(double) :: MinH, MaxH
    integer, optional :: spin ! 1 = spin up, and 2 = spin down

    ! Local variables
    integer :: np, nat, nb, ist, i, j, gcspart, wheremat, iprim, nsfi, nsfj, stat
    real(double), allocatable, dimension(:) :: HMin, HMax
    real(double) :: Hval
    logical :: do_spin_down

    do_spin_down = .false.
    if (present (spin)) then
       if (spin == 2) do_spin_down = .true. 
    end if

    MinH = 1.0e30_double
    MaxH = -1.0e30_double
    iprim = 0
    allocate(Hmin(maxnsf),Hmax(maxnsf),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating in GetHLimits: ",stat,nsfi)
    call start_timer(tmr_std_matrices)
    do np = 1, bundle%groups_on_node
       if(bundle%nm_nodgroup(np)>0) then
          do nat = 1,bundle%nm_nodgroup(np)
             iprim=iprim+1
             nsfi = mat(np,Hrange)%ndimi(nat)
             Hmin = zero
             Hmax = zero
             do nb = 1, mat(np,Hrange)%n_nab(nat)
                ist = mat(np,Hrange)%i_acc(nat) + nb - 1
                nsfj = mat(np,Hrange)%ndimj(ist)
                gcspart = BCS_parts%icover_ibeg(mat(np,Hrange)%i_part(ist)) + mat(np,Hrange)%i_seq(ist) - 1
                !wheremat = matrix_pos(matH,iprim,halo(Hrange)%i_halo(gcspart),1,1)
                wheremat = matrix_pos(matH,iprim,halo(Hrange)%i_halo(gcspart))
                if (wheremat==mat(np,Hrange)%onsite(nat)) then
                   if (do_spin_down) then
                      do i = 1, nsfj
                         do j = 1, nsfi
                            Hval = return_matrix_value (matH_dn,np,nat,iprim,nb,j,i)
                            !if(i==j) then
                            HMin(i) = HMin(i) + abs(Hval)
                            !else
                            !   HMin(i) = HMin(i) - abs (return_matrix_value (matH,np,nat,iprim,nb,j,i))
                            !end if
                            HMax(i) = HMax(i) + abs(Hval)
                         enddo
                      enddo
                   else ! if (do_spin_down) then
                      do i = 1, nsfj
                         do j = 1, nsfi
                            Hval = return_matrix_value (matH,np,nat,iprim,nb,j,i)
                            HMin(i) = HMin(i) + abs(Hval)
                            HMax(i) = HMax(i) + abs(Hval)
                         end do
                      end do
                   end if ! if (do_spin_down) then
                else ! if (wheremat==mat(np,Hrange)%onsite(nat))
                   if (do_spin_down) then
                      do i = 1,nsfj
                         do j=1,nsfi
                            Hval = return_matrix_value (matH_dn,np,nat,iprim,nb,j,i)
                            HMin(i) = HMin(i) - abs(Hval)
                            HMax(i) = HMax(i) + abs(Hval)
                         end do
                      end do
                   else ! if (do_spin_down) then
                      do i = 1,nsfj
                         do j=1,nsfi
                            Hval = return_matrix_value (matH,np,nat,iprim,nb,j,i)
                            HMin(i) = HMin(i) - abs(Hval)
                            HMax(i) = HMax(i) + abs(Hval)
                         end do
                      end do
                   end if ! if (do_spin_down) then
                endif ! if (wheremat==mat(np,Hrange)%onsite(nat))
             enddo ! do nb = 1, mat(np,Hrange)%n_nab(nat)
             ! do a simple search for the min amd max among HMin(i) and HMax(i)
             do i=1,nsfi
                if(HMin(i)<MinH) MinH = HMin(i)
                if(HMax(i)>MaxH) MaxH = HMax(i)
             enddo
          end do ! nat = bundle%nm_nodgroup
       endif
    enddo
    call stop_timer(tmr_std_matrices)
    deallocate(Hmin,Hmax)
    return
  end subroutine GetHLimits

end module McWeeny
