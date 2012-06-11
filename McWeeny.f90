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

  use global_module, only: io_lun
  use timer_stdclocks_module, only: start_timer, stop_timer, tmr_std_matrices

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
  !!   2012/03/11 L.Tong
  !!   - Major rewrite for spin implementation
  !!   - implementation of spin non-polarised case and spin fixed-case
  !!     share same piece of code. Only implementation for the variable
  !!     spin case differs, due to need to use direct sum.
  !!  SOURCE
  !!
  subroutine McWMin(n_L_iterations, deltaE)

    use datatypes
    use numbers
    use matrix_data,   only: Lrange
    use mult_module,   only: allocate_temp_matrix, matrix_sum,         &
                             free_temp_matrix, matL, matH, matS, matT, &
                             symmetrise_L,  matLS, matSL
    use global_module, only: iprint_DM, IPRINT_TIME_THRES1,            &
                             IPRINT_TIME_THRES2,                       &
                             nspin, spin_factor,                       &
                             flag_fix_spin_population
    use GenComms,      only: inode, ionode
    use timer_module,  only: cq_timer, start_timer, stop_print_timer,  &
                             WITH_LEVEL

    implicit none

    ! Passed variables
    integer,      intent(in)  :: n_L_iterations
    real(double), intent(out) :: deltaE
    ! Local variables
    integer                        :: n_iterations, spin
    integer,      dimension(nspin) :: matRhoNew, mat_oldL
    real(double), dimension(nspin) :: cn, oldE, omega1, c_old
    real(double)                   :: cn_tot, oldE_tot, omega1_tot, c_old_tot
    type(cq_timer)                 :: tmr_l_tmp1,tmr_l_tmp2
    logical,      dimension(nspin) :: done

    ! real(double) :: cn, oldE, omega1, c_old
    ! real(double) :: cn_up, cn_dn, oldE_up, oldE_dn, omega1_up, &
    !      omega1_dn, c_old_up, c_old_dn

    call start_timer(tmr_l_tmp1,WITH_LEVEL)

    done(:) = .false.
    oldE(:) = 1.0e30_double
    oldE_tot = spin_factor * sum(oldE(:))
    c_old(:) = zero
    c_old_tot = zero

    do spin = 1, nspin
       mat_oldL(spin) = allocate_temp_matrix(Lrange,0)
       matRhoNew(spin) = allocate_temp_matrix(Lrange,0)
       call matrix_sum(zero, mat_oldL(spin), one, matL(spin))
    end do

    if (nspin == 1 .or. flag_fix_spin_population) then

       iter: do n_iterations = 1, n_L_iterations

          call start_timer(tmr_l_tmp2, WITH_LEVEL)
          if (inode == ionode .and. iprint_DM >= 1) &
               write (io_lun, 1) n_iterations
          call McW_matrix_multiply_nonds(matH, matS, matL, matT, &
                                         matRhoNew, omega1, cn)

          do spin = 1, nspin
             if (omega1(spin) > oldE(spin) .or.                       &
                 abs(omega1(spin) - oldE(spin)) < 1.0e-08_double .or. &
                 cn(spin) > one .or.                                  &
                 cn(spin) < zero) then
                
                done(spin) = .true.
                
             end if
          end do
          if (done(1) .and. done(nspin)) then
             if (inode == ionode .and. iprint_DM > 0) then
                do spin = 1, nspin
                   write (io_lun, 2) spin, c_old(spin), oldE(spin),&
                                     cn(spin), omega1(spin)
                end do
             end if
             call stop_print_timer(tmr_l_tmp2, "a McWeeny iteration", &
                                   IPRINT_TIME_THRES2)
             exit iter
          end if

          deltaE = zero
          do spin = 1, nspin
              ! If new energy is lower, copy new rho into old rho and repeat
              deltaE = deltaE + spin_factor * (omega1(spin) - oldE(spin))
              call matrix_sum(zero, matL(spin), one, matRhoNew(spin))
              oldE(spin) = omega1(spin)
              c_old(spin) = cn(spin)
              if (inode == ionode .and. iprint_DM >= 0) &
                   write (io_lun, 3) spin, oldE(spin), c_old(spin)
           end do ! spin
           call symmetrise_L()
           call stop_print_timer(tmr_l_tmp2, "a McWeeny iteration",&
                                 IPRINT_TIME_THRES2)

        end do iter ! n_iterations

     else ! variable spin 

        ! If the spin populations are allowed to vary there will be
        ! only one omega1_tot and cn_tot etc for the entire direct sum

        do n_iterations = 1, n_L_iterations
           call start_timer(tmr_l_tmp2, WITH_LEVEL)
           if (inode == ionode .and. iprint_DM >= 1) &
                write(io_lun, 1) n_iterations
           ! find new rho and energy
           call McW_matrix_multiply_ds(matH, matS, matL, matT, &
                                       matRhoNew, omega1_tot, cn_tot)
           if (omega1_tot > oldE_tot .or.                      &
               abs(omega1_tot - oldE_tot) < 1.0e-8_double .or. &
               cn_tot > one .or. cn_tot < zero) then
              
              if (inode == ionode .and. iprint_DM > 0) &
                   write (io_lun, 4) c_old_tot, oldE_tot, cn_tot, omega1_tot
              call stop_print_timer(tmr_l_tmp2, "a McWeeny iteration", &
                                    IPRINT_TIME_THRES2)
              exit

           end if
           ! If new energy is lower, copy new rho into old rho and repeat
           deltaE = omega1_tot - oldE_tot
           do spin = 1, nspin
              call matrix_sum(zero, matL(spin), one, matRhoNew(spin))
           end do
           oldE_tot = omega1_tot
           c_old_tot = cn_tot
           if (inode == ionode .and. iprint_DM >= 0) &
                write (io_lun, 5) oldE_tot, c_old_tot
           call symmetrise_L()
           call stop_print_timer(tmr_l_tmp2, "a McWeeny iteration", &
                                 IPRINT_TIME_THRES2)
        end do ! n_iterations

     end if ! (nspin == 1 .or. flag_fix_spin_population)
     
     do spin = nspin, 1, -1
        if (n_iterations < 3) &
             call matrix_sum(zero, matL(spin), one, mat_oldL(spin))
        call free_temp_matrix(matRhoNew(spin))
        call free_temp_matrix(mat_oldL(spin))
     end do

     if (inode == ionode .and. iprint_DM >= 1) then
        write (io_lun, 6), n_iterations
        if (nspin == 1 .or. flag_fix_spin_population) then
           do spin = 1, nspin
              write (io_lun, 7) spin, omega1(spin)
           end do
        else
           write (io_lun, 8) omega1_tot
        end if
     end if

     call stop_print_timer(tmr_l_tmp1,"MCWEENY",IPRINT_TIME_THRES1)
     
     return

1   format(/,20x,'McWeeny L iteration:',i5)
2   format(2x,'Rounding error found (spin=',i1,'): ',4f25.15)
3   format(2x,'Energy (spin=',i1,'): ',f25.15,' Mid-point: ',f25.15)
4   format(2x,'Rounding error found (direct sum): ',4f25.15)
5   format(2x,'Energy (direct sum): ',f25.15,' Mid-point: ',f25.15)
6   format(/,20x,'Functional value reached after ',i5,' L iterations')
7   format(20x,'Omega (spin=',i1,'): ',f25.15)
8   format(20x,'Omega (direct sum): ', f25.15)

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
  !!    Bug fix - changed location of output for mubar and added gmin
  !!    and gmax calls on hmin and hmax
  !!   11:08, 04/02/2003 drb 
  !!    Implemented quick fix for bracketing mu with Hamiltonian
  !!    limits (required for stability)
  !!   12:30, 2004/06/09 dave
  !!    Fixed assumption that T and S have same range
  !!   2004/10/29 drb
  !!    Removed gsum on SX: MAJOR bug found by TM which was messing up
  !!    initialisation; also removed loops to change hmin and hmax
  !!    (unnecessary, I think)
  !!   2011/07/25 L.Tong 
  !!    Adding Spin polarisation
  !!   Friday, 2011/08/05 L.Tong 
  !!    Removed the obsolete number_of_bands parameter, and the total
  !!    electron numbers are obtained directly from the global_module
  !!   2012/01/30 L.Tong
  !!    Rewrite for the variable spin case, works for the whole direct
  !!    sum.
  !!   2012/03/08 L.Tong
  !!    Rewrote spin implementation
  !!   2012/05/27 L.Tong
  !!   - Added check for the case ne(spin) > n_o(spin), i.e. the
  !!     number of electrons requested by the user at initialisation
  !!     in one spin channel exceeds the number of orbitals avaliable
  !!     for that spin channel. The end result of this is that
  !!     eigenvalues of the L will no longer be contained within 0 and
  !!     1, due to forced double occupancy of a same spin orbital. A
  !!     warning will be produced if this happens.
  !!  SOURCE
  !!
  subroutine InitMcW(inode, ionode)

    use datatypes
    use numbers
    use global_module,  only: iprint_DM, ni_in_cell,                  &
                              IPRINT_TIME_THRES1, nspin,              &
                              flag_fix_spin_population, ne_in_cell,   &
                              ne_spin_in_cell, spin_factor
    use species_module, only: species, nsf_species
    use mult_module,    only: allocate_temp_matrix, matH, matS, matT, &
                              matL, mult, matrix_trace, T_S_TS,       &
                              matrix_product_trace, matrix_product,   &
                              free_temp_matrix
    use matrix_data,    only: Lrange, Srange, TSrange
    use GenBlas
    use GenComms,       only: gsum, my_barrier, gmin, gmax
    use timer_module,   only: cq_timer,start_timer, stop_print_timer, &
                              WITH_LEVEL

    implicit none

    ! Passed variables
    ! real(double) :: number_of_bands 
    integer :: inode, ionode

    ! Local variables
    real(double) :: n_e, n_o, hmax_ds, hmin_ds, SX, SXHX_ds, A
    real(double) :: direct_sum_factor = one

    real(double), dimension(nspin) :: ne, mu1, mu2, mubar, lambda, &
                                      hmin, hmax, SXHX
    integer,      dimension(nspin) :: matXHX, matSXHX, mat_temp
    integer                        :: matTS, length, i, spin
    type(cq_timer)                 :: tmr_l_tmp1

    call start_timer(tmr_l_tmp1, WITH_LEVEL)
    if (inode == ionode .and. iprint_DM > 1) write (io_lun, 1)
    ! We must first initialise rho
    ! n_e is the number of electrons per spin channel
    ne(1:nspin) = ne_spin_in_cell(1:nspin)
    n_e = spin_factor * sum(ne(:))
    ! n_e is total number of electrons

    ! we are working with the whole direct sum for variable spin case
    ! so we need to have twice the size of basis set
    if (nspin == 2 .and. (.not. flag_fix_spin_population)) &
         direct_sum_factor = two

    ! n_o is number of support functions 
    n_o = zero
    do i = 1, ni_in_cell
       n_o = n_o + direct_sum_factor * nsf_species(species(i))
    end do
    
    ! check if we have enough orbitals for the given electron
    ! populations for the initialisation procedure to work properly.
    do spin = 1, nspin
       if (ne(spin) / n_o * direct_sum_factor > one) then
          write (io_lun, '(/,2x,a,i1,a,i1,a,f12.6,a,f12.6,a,/)')      &
               '*** WARNING!!! You have too few number of support &
                &functions for the calculation. The number of &
                &electrons in spin = ', spin, ' is greater than the &
                &number of orbitals!. ne(', spin, ') = ', ne(spin),   &
                ', n_o for spin channel = ', n_o / direct_sum_factor, &
                ' ***'
       end if
    end do

    if (nspin == 2) then
       if (inode == ionode .and. iprint_DM > 1) &
            write (io_lun, 3) ne(1), ne(2), n_o
    else
       if (inode == ionode .and. iprint_DM > 1) &
            write (io_lun, 2) n_e, n_o
    end if

    matTS = allocate_temp_matrix(TSrange,0)
    call matrix_product(matT, matS, matTS, mult(T_S_TS))
    SX = direct_sum_factor * matrix_trace(matTS)
    if (inode == ionode .and. iprint_DM > 1) then 
       write (io_lun, fmt='(2x,"SX is    ",f25.15)') SX
       write (io_lun, fmt='(2x,"SX/no is ",f25.15)') SX / n_o
    end if

    SXHX_ds = zero
    do spin = 1, nspin
       matXHX(spin) = allocate_temp_matrix(Lrange, 0)
       mat_temp(spin) = allocate_temp_matrix(Srange, 0)
       call McWXHX(matH(spin), matT, matXHX(spin), mat_temp(spin))
       SXHX(spin) = matrix_product_trace(matS, mat_temp(spin))
       SXHX_ds = SXHX_ds + spin_factor * SXHX(spin)
       if (inode == ionode .and. iprint_DM > 1) &
            write (io_lun, '(2x,"SXHX(spin=",i1,") = ",f25.15)') &
                  spin, SXHX(spin)
    end do
    if (inode == ionode .and. iprint_DM > 1) &
            write (io_lun, '(2x,"SXHX_ds = ",f25.15)') &
                  SXHX_ds    

    call my_barrier()

    ! Get Gershgorin limits on H
    call GetHLimits(hmin, hmax)
    do spin = 1, nspin
       call gmin(hmin(spin))
       call gmax(hmax(spin))
       if (inode == ionode .and. iprint_DM > 1) &
          write (io_lun, 4) spin, hmin(spin), hmax(spin)  
    end do
    ! for variable spin case find min and max of the direct sum
    if (nspin == 2 .and. (.not. flag_fix_spin_population)) then
       hmin_ds = minval(hmin(:))
       hmax_ds = maxval(hmax(:))
       if (inode == ionode .and. iprint_DM > 1) write (io_lun, 5) hmin, hmax
    end if
    
    A = (one - SX/n_o)

    if (nspin == 1 .or. flag_fix_spin_population) then
       do spin = 1, nspin
          mu1(spin) = SXHX(spin) / n_o + hmax(spin) * A
          mu2(spin) = (hmin(spin) * A * ne(spin) / (n_o - ne(spin)) - &
                       SXHX(spin) / n_o) / &
                      (n_o * A / (n_o - ne(spin)) - one)
          if (inode == ionode) &
               write (io_lun, &
                      '(2x,"Mu1, Mu2, for spin = ",i1," are: ",2f25.15)') &
                     spin, mu1(spin), mu2(spin)
          if ((mu1(spin) < hmax(spin) .and. mu1(spin) > hmin(spin)) .and. &
              (abs(ne(spin) / (hmax(spin) - mu1(spin))) < &
               abs((n_o - ne(spin)) / (mu2(spin) - hmin(spin))))) then
             mubar(spin) = mu1(spin)
             lambda(spin) = ne(spin) / (hmax(spin) - mu1(spin))
          else
             mubar(spin) = mu2(spin)
             lambda(spin) = (n_o - ne(spin)) / (mu2(spin) - hmin(spin))
          endif
          if (inode == ionode .and. iprint_DM > 1) &
               write (io_lun, &
                      '(2x,"mubar, lambda for spin = ",i1," are: ",2f25.15)') &
                     spin, mubar(spin), lambda(spin)
       end do
       ! Calculate L0. 
       do spin = 1, nspin
          call McWRho0(matL(spin), matT, matXHX(spin), mubar(spin), &
               lambda(spin), ne(spin), n_o)
       end do
    else ! variable spin
       ! mu1(1) and mu1(2) always equal in this case. Same for mu2
       mu1(:) = SXHX_ds / n_o + hmax_ds * A
       mu2(:) = (hmin_ds * A * n_e / (n_o - n_e) - SXHX_ds / n_o) / &
                (n_o * A / (n_o - n_e) - one)
       if (inode == ionode) &
            write (io_lun, '(2x,"Mu1, Mu2: ",2f25.15)') mu1(1), mu2(1)
       if ((mu1(1) < hmax_ds .and. mu1(1) > hmin_ds) .and. &
            (abs(n_e / (hmax_ds - mu1(1))) < &
             abs((n_o - n_e) / (mu2(1) - hmin_ds)))) then
          mubar(:) = mu1(:)
          lambda(:) = n_e / (hmax_ds - mu1(:))
       else
          mubar(:) = mu2(:)
          lambda(:) = (n_o - n_e) / (mu2(:) - hmin_ds)
       endif
       if (inode == ionode .and. iprint_DM > 1) then
          write (io_lun, 6) mubar(1)
          write (io_lun, '(2x,"lambda is ",f25.15)') lambda(1)
       end if
       ! Calculate L0. 
       do spin = 1, nspin
          call McWRho0(matL(spin), matT, matXHX(spin), mubar(spin), &
                       lambda(spin), n_e, n_o)
       end do
    end if
    
    ! Free the temporary matrices, order of deallocation is important
    do spin = nspin, 1, -1
       call free_temp_matrix(mat_temp(spin))
       call free_temp_matrix(matXHX(spin))
    end do
    call free_temp_matrix(matTS)

    call stop_print_timer(tmr_l_tmp1, "McWeeny initialisation", &
                          IPRINT_TIME_THRES1)

1   format(1x,'Welcome to InitMcW')
2   format(2x,'Electrons: ',f15.6,' Orbitals: ',f15.6)
3   format(2x,'Electrons_up: ',f15.6,' Electrons_dn: ',f15.6,' Orbitals: ',f15.6)
4   format(2x,'Minimum and maximum limites on H(spin=',i1') are ',2f15.6)
5   format(2x,'Minimum and maximum limites on overall H are ',2f15.6)
6   format(2x,'Mubar is ',f15.6)

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
  !!   2012/03/11 L.Tong
  !!    - Major rewrite of spin implementation
  !!    - Use matLS, matSL from mult module directly
  !!    - All inputs except matS and matInvS are arrays of dimension nspin
  !!  SOURCE
  !!
  subroutine McW_matrix_multiply_nonds(matH, matS, matL, matInvS, &
                                       matRhoNew, energy, c)
    use datatypes
    use numbers
    use primary_module, only: bundle
    use matrix_data,    only: LSLrange, Lrange, Hrange, Srange
    use mult_module,    only: matrix_product, matrix_sum,             &
                              matrix_product_trace, matrix_scale,     &
                              matrix_transpose, allocate_temp_matrix, &
                              free_temp_matrix, L_S_LS, LS_L_LSL,     &
                              LSL_SL_L, mult, matLS, matSL
    use global_module,  only: iprint_DM, nspin, spin_factor
    use GenComms,       only: gsum, cq_abort, inode, ionode, my_barrier

    implicit none

    ! Passed variables
    integer,      dimension(:) :: matH, matL, matRhoNew
    integer                    :: matS, matInvS
    real(double), dimension(:) :: energy, c
    ! Local Variables
    integer :: spin
    real(double), dimension(nspin) :: c1, c2, cn
    integer,      dimension(nspin) :: matLSL, matLSLSL, mat_top, mat_bottom

    do spin = 1, nspin
       call matrix_product(matL(spin), matS, matLS(spin), mult(L_S_LS))
       ! Generate LSL, LSLSL
       matLSL(spin) = allocate_temp_matrix(LSLrange,0)
       matLSLSL(spin) = allocate_temp_matrix(Lrange,0)
       mat_top(spin) = allocate_temp_matrix(Srange,0)
       mat_bottom(spin) = allocate_temp_matrix(Srange,0)
       call my_barrier
       call matrix_product(matLS(spin), matL(spin), matLSL(spin), &
                           mult(LS_L_LSL))
       call matrix_product(matLSL(spin), matLS(spin), matLSLSL(spin), &
                           mult(LSL_SL_L))
       !----------------------------------------------------------------
       ! Shorten L, LSL and LSLSL to range S and take traces (can do it
       ! with sdot because S is symmetric)
       !----------------------------------------------------------------
       call matrix_sum(zero, mat_top(spin),  one, matLSL(spin))
       call matrix_sum(one,  mat_top(spin), -one, matLSLSL(spin))
       c1(spin) = matrix_product_trace(matS, mat_top(spin))
       if (inode == ionode .and. iprint_DM >= 2) &
            write (io_lun, '(2x,"S.top (spin=",i1,") is  ",f25.15)') &
                  spin, c1(spin)
       call matrix_sum(zero, mat_bottom(spin), one, matL(spin))
       c1(spin) = matrix_product_trace(matS, mat_bottom(spin))
       if (inode == ionode .and. iprint_DM >= 2) &
            write (io_lun, '(2x,"N_e (spin=",i1,") is  ",f25.15)') &
                  spin, c1(spin)
       call matrix_sum(one, mat_bottom(spin), -one, matLSL(spin))
       c1(spin) = matrix_product_trace(matS, mat_top(spin))
       c2(spin) = matrix_product_trace(matS, mat_bottom(spin))
       if (c2(spin) /= zero) then
          cn(spin) = c1(spin) / c2(spin)
       else
          call cq_abort('McW_matrix_multiply: c2 is zero for spin', spin)
       endif
       c(spin) = cn(spin)
       if (inode == ionode .and. iprint_DM >= 2) &
            write (io_lun, '(2x,"c, c1, c2 (spin=",i1,") are ",3f25.15)') &
                  spin, cn(spin), c1(spin), c2(spin)                              
       ! Shorten LSL and LSLSL to range L
       call matrix_scale(zero, matRhoNew(spin))
       call matrix_sum(zero, matRhoNew(spin), -one, matLSLSL(spin))
       call matrix_sum(one, matRhoNew(spin), one + cn(spin), matLSL(spin))
       if(cn(spin) < half) then
          call matrix_sum(one, matRhoNew(spin), (one - two * &
                          cn(spin)), matL(spin))
          call matrix_scale((one / (one - cn(spin))), matRhoNew(spin))
       else
          call matrix_scale(one / cn(spin), matRhoNew(spin))
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
       !%%! if(inode.eq.ionode.AND.iprint_DM>=2) &
       !%%!     write(io_lun,*) 'Idempot monitor is ',cn/256.0_double

       ! Shorten rhonew to range H into temp
       call free_temp_matrix(mat_bottom(spin))
       mat_bottom(spin) = allocate_temp_matrix(Hrange, 0)
       call matrix_sum(zero, mat_bottom(spin), one, matRhoNew(spin))
       energy(spin) = matrix_product_trace(mat_bottom(spin), matH(spin))
       if (inode == ionode .and. iprint_DM >= 3) &
            write (io_lun, '(2x,"energy (spin=",i1,") is ",f25.15)') &
                  spin, energy(spin)
       call matrix_sum(zero, mat_top(spin), one, matRhoNew(spin))
       c1(spin) = matrix_product_trace(mat_top(spin), matS)
       if (inode == ionode .and. iprint_DM >= 3) &
            write (io_lun, '(2x,"N_e(2) (spin=",i1,") is ",f25.15)') &
                  spin, c1(spin)
       ! note that N_e is printed without factor of two here even for
       ! spin un-polarised results.
       call free_temp_matrix(mat_bottom(spin))
       call free_temp_matrix(mat_top(spin))
       call free_temp_matrix(matLSLSL(spin))
       call free_temp_matrix(matLSL(spin))
       if (abs(c2(spin)) < 1.0e-6) c(spin) = -0.01_double
    end do ! spin

    return

  end subroutine McW_matrix_multiply_nonds
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
  !!   2012/03/11 L.Tong
  !!   - Major rewrite of the spin implementation
  !!   - all inputs except matS and matInvS are arrays of dimension nspin
  !! SOURCE
  !!
  subroutine McW_matrix_multiply_ds(matH, matS, matL, matInvS, &
                                    matRhoNew, energy_total, c_total)

    use datatypes
    use numbers
    use primary_module, only: bundle
    use matrix_data,    only: LSLrange, Lrange, Hrange, Srange
    use mult_module,    only: matrix_product, matrix_sum,             &
                              matrix_product_trace, matrix_scale,     &
                              matrix_transpose, allocate_temp_matrix, &
                              free_temp_matrix, matLS, matSL, L_S_LS, &
                              LS_L_LSL, LSL_SL_L, mult
    use global_module,  only: iprint_DM, nspin, spin_factor
    use GenComms,       only: gsum, cq_abort, inode, ionode, my_barrier

    implicit none

    ! Passed variables
    integer, dimension(:) :: matH, matL, matRhoNew
    integer               :: matS, matInvS
    real(double)          :: energy_total, c_total

    ! Local Variables
    real(double)              :: c1, c2, tmp
    integer, dimension(nspin) :: matLSL, matLSLSL, mat_top, mat_bottom
    integer                   :: spin
    
    c1 = zero
    c2 = zero
    energy_total = zero

    do spin = 1, nspin
       call matrix_product(matL(spin), matS, matLS(spin), mult(L_S_LS))
       ! Generate LSL, LSLSL
       matLSL(spin) = allocate_temp_matrix(LSLrange, 0)
       matLSLSL(spin) = allocate_temp_matrix(Lrange, 0)
       mat_top(spin) = allocate_temp_matrix(Srange, 0)
       mat_bottom(spin) = allocate_temp_matrix(Srange, 0)
       call my_barrier
       call matrix_product(matLS(spin), matL(spin), matLSL(spin), &
                           mult(LS_L_LSL))
       call matrix_product(matLSL(spin), matLS(spin), matLSLSL(spin), &
                           mult(LSL_SL_L))
       !----------------------------------------------------------------
       ! Shorten L, LSL and LSLSL to range S and take traces (can do it
       ! with sdot because S is symmetric)
       !----------------------------------------------------------------
       call matrix_sum(zero, mat_top(spin),  one, matLSL(spin))
       call matrix_sum(one,  mat_top(spin), -one, matLSLSL(spin))
       tmp = matrix_product_trace(matS, mat_top(spin))
       if (inode == ionode .and. iprint_DM >= 2) &
            write (io_lun, '(2x,"S.top (spin=",i1,") is ",f25.15)') spin, tmp
       call matrix_sum(zero, mat_bottom(spin), one, matL(spin))
       tmp = matrix_product_trace(matS, mat_bottom(spin))
       if (inode == ionode .and. iprint_DM >= 2) &
            write (io_lun, '(2x,"N_e (spin=",i1,") is ",f25.15)') spin, tmp
       ! get the real mat_bottom for spin up
       call matrix_sum(one, mat_bottom(spin), -one, matLSL(spin))
       ! c1 and c2 calculated
       c1 = c1 + spin_factor * matrix_product_trace(matS, mat_top(spin))
       c2 = c2 + spin_factor * matrix_product_trace(matS, mat_bottom(spin))
       ! finished  mat_bottom here, change to Hrange for something else
       call free_temp_matrix(mat_bottom(spin))
       mat_bottom(spin) = allocate_temp_matrix(Hrange,0)
    end do ! spin
    if (c2 /= zero) then
       c_total = c1 / c2
    else
       call cq_abort('McW_matrix_multiply_ds: c2 is zero')
    end if
    if (inode == ionode .and. iprint_DM >= 2) &
         write (io_lun, '(2x,"c_total, c1, c2 are ",3f25.15)') &
               c_total, c1, c2
    ! Shorten LSL and LSLSL to range L
    do spin = 1, nspin
       call matrix_scale(zero, matRhoNew(spin))
       call matrix_sum(zero, matRhoNew(spin), -one, matLSLSL(spin))
       call matrix_sum(one, matRhoNew(spin), one + c_total, matLSL(spin))
       if (c_total < half) then
          call matrix_sum(one, matRhoNew(spin), (one - two * c_total),&
                          matL(spin))
          call matrix_scale((one/(one - c_total)), matRhoNew(spin))
       else
          call matrix_scale(one/c_total, matRhoNew(spin))
       end if
       call matrix_sum(zero, mat_bottom(spin), one, matRhoNew(spin))
       energy_total = energy_total + spin_factor * &
                      matrix_product_trace(mat_bottom(spin), matH(spin))
    end do ! spin
    if (inode == ionode .and. iprint_DM >= 3) &
         write (io_lun, '(2x,"energy_total: ",f25.15)') energy_total
    do spin = 1, nspin
       call matrix_sum(zero, mat_top(spin), one, matRhoNew(spin))
       tmp = matrix_product_trace(mat_top(spin), matS)
       if (inode == ionode .and. iprint_DM >= 3) &
            write (io_lun, '(2x,"N_e (spin=",i1,") is ",f25.15)') &
                  spin, tmp
    end do ! spin

    ! free temp matrices
    do spin = nspin, 1, -1
       call free_temp_matrix(mat_bottom(spin))
       call free_temp_matrix(mat_top(spin))
       call free_temp_matrix(matLSLSL(spin))
       call free_temp_matrix(matLSL(spin))
    end do

    return

  End Subroutine McW_matrix_multiply_ds
  !!***


  !!****f* MacWeeny/McWXHX
  !! PURPOSE
  !!  sbrt McWXHX: Pre- and post- multiplies H by S^-1, and returns
  !!  the result to both L range (bab) and S range (tm)
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
    use matrix_data,   only: THrange
    use mult_module,   only: matrix_product, matrix_sum, T_H_TH, &
                             TH_T_L, allocate_temp_matrix, mult, &
                             free_temp_matrix, matrix_trace
    use GenComms,      only: inode, ionode
    use global_module, only: iprint_DM

    implicit none

    ! Passed variables
    integer :: matA, matB, matBAB, mat_temp

    ! Local variable
    integer :: matBA
    real(double) :: xx

    matBA = allocate_temp_matrix(THrange,0)
    call matrix_product(matB, matA, matBA, mult(T_H_TH))
    if (iprint_DM >= 2) xx = matrix_trace(matBA)
    if (inode == ionode .and. iprint_DM >= 2) &
         write (io_lun, '(2x,"Trace of BA:            ",f25.15)') xx
    call matrix_product(matBA, matB, matBAB, mult(TH_T_L))
    if (iprint_DM >= 2) xx = matrix_trace(matBAB)
    if (inode == ionode .and. iprint_DM >= 2) &
         write (io_lun, '(2x,"Trace of BAB (L range): ",f25.15)') xx
    call matrix_sum(zero, mat_temp, one, matBAB)
    if (iprint_DM >= 2) xx = matrix_trace(mat_temp)
    if (inode == ionode .and. iprint_DM >= 2) &
         write (io_lun, '(2x,"Trace of BAB (S range): ",f25.15)') xx
    call free_temp_matrix(matBA)
    return
  end Subroutine McWXHX
  !!***


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
    use GenComms,      only: inode, ionode
    use mult_module,   only: matrix_sum
    use global_module, only: iprint_DM

    implicit none

    integer :: matA, matB, matC
    Real( double ) :: m, l, n_e, n_o, tmp

    tmp = (l*m+n_e)/n_o
    if (inode == ionode .and. iprint_DM >= 2) &
         write (io_lun, '(2x,"l, m, n_e, n_o and tmp are ",5f25.15)') &
               l, m, n_e, n_o, tmp
    call matrix_sum(zero, matA, tmp, matB)
    tmp = -l/n_o
    if (inode == ionode .and. iprint_DM >= 2) &
         write (io_lun, '(2x,"tmp is ",f25.15)') tmp
    call matrix_sum(one, matA, tmp, matC)
    return
  end Subroutine McWRho0
  !!***

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
  !!   2012/03/11 L.Tong
  !!   - Major rewrite of spin implementation
  !!   - removed optional parameter spin
  !!   - MinH and MaxH are now arrays of spin, each contain min and max
  !!     of H in each spin channels.
  !! SOURCE
  !!
  subroutine GetHLimits(MinH, MaxH)

    use datatypes
    use numbers
    use matrix_module,  only: matrix, matrix_halo
    use matrix_data,    only: mat, Hrange, halo
    use cover_module,   only: BCS_parts
    use primary_module, only: bundle
    use mult_module,    only: matH, return_matrix_value, matrix_pos
    use GenComms,       only: cq_abort
    use maxima_module,  only: maxnsf
    use global_module,  only: nspin

    implicit none

    ! Passed variables
    real(double), dimension(:) :: MinH, MaxH

    ! Local variables
    integer :: np, nat, nb, ist, i, j, gcspart, wheremat, iprim, nsfi,&
               nsfj, stat, spin
    real(double), dimension(maxnsf,nspin) :: HMin, HMax
    real(double) :: Hval

    MinH(:) = 1.0e30_double
    MaxH(:) = -1.0e30_double

    call start_timer(tmr_std_matrices)
    do spin = 1, nspin
       iprim = 0
       do np = 1, bundle%groups_on_node
          if (bundle%nm_nodgroup(np) > 0) then
             do nat = 1,bundle%nm_nodgroup(np)
                iprim=iprim+1
                nsfi = mat(np,Hrange)%ndimi(nat)
                HMin(:,spin) = zero
                HMax(:,spin) = zero
                do nb = 1, mat(np,Hrange)%n_nab(nat)
                   ist = mat(np,Hrange)%i_acc(nat) + nb - 1
                   nsfj = mat(np,Hrange)%ndimj(ist)
                   gcspart = &
                        BCS_parts%icover_ibeg(mat(np,Hrange)%i_part(ist)) +&
                        mat(np,Hrange)%i_seq(ist) - 1
                   wheremat = matrix_pos(matH(spin), iprim, &
                                         halo(Hrange)%i_halo(gcspart))
                   if (wheremat == mat(np, Hrange)%onsite(nat)) then
                      do i = 1, nsfj
                         do j = 1, nsfi
                            Hval = return_matrix_value(matH(spin), np, &
                                                       nat, iprim, nb, j, i)
                            HMin(i,spin) = HMin(i,spin) + abs(Hval)
                            HMax(i,spin) = HMax(i,spin) + abs(Hval)
                         end do
                      end do
                   else
                      do i = 1,nsfj
                         do j= 1, nsfi
                            Hval = return_matrix_value(matH(spin), np, &
                                                       nat, iprim, nb, j, i)
                            HMin(i,spin) = HMin(i,spin) - abs(Hval)
                            HMax(i,spin) = HMax(i,spin) + abs(Hval)
                         end do
                      end do
                   end if
                end do ! do nb = 1, mat(np,Hrange)%n_nab(nat)
                ! do a simple search for the min amd max among HMin(i) and HMax(i)
                do i = 1, nsfi
                   if (HMin(i,spin) < MinH(spin)) MinH(spin) = HMin(i,spin)
                   if (HMax(i,spin) > MaxH(spin)) MaxH(spin) = HMax(i,spin)
                enddo
             end do ! nat = bundle%nm_nodgroup
          end if ! (bundle%nm_nodgroup(np) > 0)
       end do ! np
    end do ! spin
    call stop_timer(tmr_std_matrices)

    return
  end subroutine GetHLimits
  !!***

end module McWeeny
