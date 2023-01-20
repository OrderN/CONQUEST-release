! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! Module blip_minimisation
! ------------------------------------------------------------------------------
! Code area 11: basis functions
! ------------------------------------------------------------------------------

!!****h* Conquest/blip_minimisation *
!!  NAME
!!   blip_minimisation
!!  PURPOSE
!!   Hold the different routines associated with blip minimisation
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   09:17, 2003/03/24 dave
!!  MODIFICATION HISTORY
!!   10:09, 13/02/2006 drb 
!!    Removed all explicit references to data_ variables and rewrote in terms of new 
!!    matrix routines
!!   2008/02/04 08:27 dave
!!    Changed for output to file not stdout
!!   2008/06/10 ast
!!    Added timers
!!   2012/04/02 17:21 dave
!!    Changes for analytic blip integration
!!   2014/09/15 18:30 lat
!!    fixed call start/stop_timer to timer_module (not timer_stdlocks_module !)
!!  SOURCE
!!
module blip_minimisation

  use datatypes
  use global_module,          only: io_lun, area_basis
  use timer_module,           only: start_timer, stop_timer
  use timer_stdclocks_module, only: tmr_std_basis, tmr_std_allocation


  implicit none

  real(double), save :: dE_blip
!!***

contains

  !!****f* blip_minimisation/vary_support *
  !!
  !!  NAME 
  !!   vary_support
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Performs the minimisation with respect to the 
  !!    support functions. The method used is CG. This
  !!    subroutine follows closely the strategy of the
  !!    original konquest program. (Wow)
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !!   blip, common, datatypes, DiagModule, GenBlas, GenComms, logicals,
  !!   matrix_data, maxima_module, mult_module, numbers, PosTan
  !!  AUTHOR
  !!   D.R.Bowler/E.H.Hernandez
  !!  CREATION DATE
  !!   03/04/95
  !!  MODIFICATION HISTORY
  !!   updated to make use of blip functions.
  !!   updated to include electron gradient
  !!   25/1/96 CMG/EHE line minimisation moved out
  !!   19/5/97 DRB added HeadGordon
  !!   23/05/2001 dave
  !!    Shortened calls to get_blip_gradient and get_electron_gradient
  !!    and added ROBODoc header
  !!   23/05/2001 dave
  !!    Shortened call to line_minimise_support
  !!   08/06/2001 dave
  !!    Added RCS Id and Log tags and changed to use GenComms
  !!   17/06/2002 dave
  !!    Tidied headers, added check for solution method
  !!   09:18, 2003/03/24 dave
  !!    Included in blip_minimisation
  !!   08:29, 2003/04/03 dave
  !!    Changed to use blip_gradient for get_blip_gradient and
  !!    get_electron_gradient
  !!   11:30, 12/11/2004 dave 
  !!    Changed to use mx_at_prim from maxima
  !!   2007/03/29 08:20 dave and tsuyoshi
  !!   2007/04/17 09:36 dave
  !!    Added n_L_iterations
  !!   2008/06/10 ast
  !!    Added timers
  !!   2011/11/17 10:21 dave
  !!    Changes to blip data
  !!   2011/11/29 L.Tong
  !!    Added spin polarisation
  !!   2011/12/05 L.Tong
  !!    - Changed the temp array sum to summ to avoid potential confusion
  !!      with the intrinsic function sum
  !!    - Removed redundant parameter number_of_bands
  !!   2012/03/01 L.Tong
  !!    Fixed a bug associated to pulay forces. It is solved by putting
  !!    an loop around the calculation of the gradient at the end of the
  !!    support variation loop, so that if we have eached the end of the
  !!    blip minimisation cycle (without convergence) the gradient and
  !!    density don't get called for the one last time (which is
  !!    intended to prepare for the next iteration).
  !!   2012/03/23 L.Tong
  !!   - Changed spin implementation
  !!   - Removed redundant input parameter real(double) mu
  !!   2012/04/02 17:22 dave
  !!    Update for analytic blip integration
  !!   2016/07/13 18:30 nakata
  !!    Renamed H_on_supportfns -> H_on_atomfns
  !!   2016/07/29 18:30 nakata
  !!    Renamed supports_on_atom -> blips_on_atom
  !!   2016/08/08 15:30 nakata
  !!    Renamed supportfns -> atomfns
  !!   2017/02/23 dave
  !!    - Changing location of diagon flag from DiagModule to global and name to flag_diagonalisation
  !!   2022/12/14 13:03 dave
  !!    Change to save initial energy in total_energy_0 from new_SC_potl
  !!  SOURCE
  !!
  subroutine vary_support(n_support_iterations, fixed_potential, &
                          vary_mu, n_L_iterations, L_tolerance,  &
                          sc_tolerance, energy_tolerance,        &
                          total_energy_last, expected_reduction)

    use datatypes
    use logicals
    use numbers
    use units
    use mult_module,         only: LNV_matrix_multiply, matM12, matM4
    use GenBlas,             only: dot, copy, axpy
    use PosTan,              only: PulayC, PulayBeta, SCC, SCBeta
    use blip,                only: PreCond, blip_info
    use GenComms,            only: my_barrier, gsum, inode, ionode, cq_abort
    !use DiagModule,          only: diagon
    use blip_gradient,       only: get_blip_gradient, get_electron_gradient
    use global_module,       only: flag_precondition_blips, iprint_basis, &
                                   area_basis, nspin, spin_factor,        &
                                   flag_analytic_blip_int, flag_diagonalisation
    use io_module,           only: dump_blip_coeffs
    use functions_on_grid,   only: atomfns, H_on_atomfns,                 &
                                   gridfunctions, fn_on_grid
    use support_spec_format, only: coefficient_array,                     &
                                   coeff_array_size, grad_coeff_array,    &
                                   elec_grad_coeff_array,                 &
                                   blips_on_atom
    use primary_module,      only: bundle
    use memory_module,       only: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use S_matrix_module,     only: get_S_matrix
    use SelfCon,           only: new_SC_potl

    implicit none

    !     Shared variables
    logical      :: vary_mu, fixed_potential, convergence_flag
    integer      :: n_L_iterations
    integer      :: n_support_iterations
    real(double) :: expected_reduction, total_energy_last, &
                    energy_tolerance, L_tolerance, sc_tolerance

    ! Local variables
    real(double) :: tolerance, con_tolerance, dgg, gamma, gg, &
                    electrons_tot, sum_0, diff, total_energy_0, &
                    energy_in_tot, last_step, dN_dot_de, dN_dot_dN
    integer      :: length, n_iterations, n_tries, offset
    integer      :: k, i, j, n, spec, stat, spin
    logical      :: notredone, reduced, reset_L
    real(double), parameter :: gamma_max = 6.0_double !! TM 2007.03.29
    real(double), dimension(nspin)          :: electrons, energy_in
    real(double), dimension(:), allocatable :: search_direction, last_sd, Psd
    real(double), dimension(:), allocatable :: summ

    call start_timer(tmr_std_basis)

    allocate(search_direction(coeff_array_size), &
             last_sd(coeff_array_size), Psd(coeff_array_size), STAT=stat)
    if (stat /= 0) &
         call cq_abort("vary_support: Error alloc mem: ", coeff_array_size)
    call reg_alloc_mem(area_basis, 3*coeff_array_size, type_dbl)

    if (inode == ionode) &
         write (io_lun, fmt='(/6x,"Starting blip-coefficient variation",/)')
    if (inode == ionode .and. iprint_basis > 0) &
         write (io_lun, &
                fmt='(6x,"Performing at most ",i4," iterations &
                      &with a tolerance of ",e12.5,/)') &
               n_support_iterations, energy_tolerance

    ! Set tolerances for self-consistency and L minimisation
    con_tolerance = SCC * expected_reduction**SCBeta
    tolerance = PulayC * (0.1_double * expected_reduction)**PulayBeta
    if (con_tolerance < sc_tolerance) con_tolerance = sc_tolerance
    if (con_tolerance < ten * tolerance) &
         tolerance = 0.1_double * con_tolerance
    con_tolerance = sc_tolerance
    tolerance = L_tolerance

    if (inode == ionode .and. iprint_basis > 1) &
         write (io_lun, fmt='(6x,"L, SC tolerances: ",2f20.12)') &
               con_tolerance, tolerance

    length = coeff_array_size
    search_direction = zero
    Psd = zero
    grad_coeff_array = zero
    elec_grad_coeff_array = zero

    ! evaluate the action of H on the support functions, calculate the
    ! corresponding matrix elements in the 'support representation'
    ! and calculate the start-point energy

    total_energy_0 = total_energy_last
    if (total_energy_last == zero) &
         total_energy_0 = expected_reduction
    call my_barrier()
    total_energy_last = total_energy_0

    reset_L = .true.
    call new_SC_potl(.false., con_tolerance, reset_L,             &
         fixed_potential, vary_mu, n_L_iterations, &
         tolerance, total_energy_last)
    total_energy_0 = total_energy_last
    ! now obtain the gradient of the energy with respect to support
    ! functions. For diagonalisation, this is already stored
    if (inode == ionode .and. iprint_basis > 2) &
         write (io_lun, fmt='(6x, "Finding gradient")')
    if (.not. flag_diagonalisation) then
       call LNV_matrix_multiply(electrons, energy_in, doK, doM1, doM2,&
                                dontM3, doM4, dontphi, dontE,         &
                                mat_M12=matM12, mat_M4=matM4)
    end if
    if(flag_analytic_blip_int) call get_S_matrix(inode, ionode)
    call get_blip_gradient(inode, ionode)
    if (inode == ionode .and. iprint_basis > 2) &
         write (io_lun, fmt='(6x,"Got blip gradient")')
    call my_barrier()
    ! The density matrix is effectively idempotent during
    ! diagonalisation
    if (flag_diagonalisation) then
       gridfunctions(H_on_atomfns(1))%griddata = zero
    else
       call get_electron_gradient(atomfns, H_on_atomfns(1), &
                                  inode, ionode)
       call my_barrier()
       if (inode == ionode .AND. iprint_basis > 2) &
            write (*,fmt='(6x,"Got electron gradient")')
    end if

    ! Now we have a basic gradient, so loop
    dgg = zero
    last_step = 1.0e10_double

    ! now loop over search directions
    do n_iterations = 1, n_support_iterations

       if (inode == ionode) then
          if (iprint_basis > 0) then
             write (io_lun, 7) n_iterations
          else 
             write (io_lun, fmt='(6x,"Support variation: ",i5)') &
                   n_iterations
          end if
       end if
       ! We need the last search direction for CG manipulations
       call copy(length, search_direction, 1, last_sd, 1)

       ! The basis for searching is data_dblip
       call copy(length, grad_coeff_array, 1, search_direction, 1)
       ! Now project dblip tangential to the constant Ne hyperplane
       if (.not. flag_diagonalisation) then
          dN_dot_de = dot(length, grad_coeff_array, 1, &
                          elec_grad_coeff_array, 1)
          dN_dot_dN = dot(length, elec_grad_coeff_array, 1, &
                          elec_grad_coeff_array, 1)
          call gsum(dN_dot_de)
          call gsum(dN_dot_dN)
          if (inode == ionode .and. iprint_basis > 2) &
               write (io_lun, fmt='(6x,"dN.de, dN.dN ",2f20.12)') &
                     dN_dot_de, dN_dot_dN
          call axpy(length, - (dN_dot_de / dN_dot_dN), &
                    elec_grad_coeff_array, 1, search_direction, 1 )      
       end if
       ! Precondition if we need to
       if (flag_precondition_blips) then
          offset = 0
          if (inode == ionode .and. iprint_basis > 3) &
               write (io_lun, fmt='(6x,"Preconditioning !")')
          call my_barrier()
          do i = 1, bundle%n_prim
             spec = bundle%species(i)
             call start_timer(tmr_std_allocation)
             allocate(summ(blips_on_atom(i)%nsuppfuncs))
             call stop_timer(tmr_std_allocation)
             do j = 1, blip_info(spec)%NBlipsRegion ! In future, base on species
                summ = zero
                do k = 1, blip_info(spec)%NBlipsRegion ! Again, base on species
                   do n = 1, blips_on_atom(i)%nsuppfuncs
                      summ(n) = summ(n) + PreCond(spec)%coeffs(j,k) * &
                                search_direction(offset + (n-1) *     &
                                blip_info(spec)%NBlipsRegion + k)
                   end do
                end do
                do n = 1, blips_on_atom(i)%nsuppfuncs
                   Psd(offset + (n - 1) * blip_info(spec)%NBlipsRegion + j) = &
                        summ(n)
                end do
             end do
             offset = offset + blips_on_atom(i)%nsuppfuncs * &
                      blip_info(spec)%NBlipsRegion
             call start_timer(tmr_std_allocation)
             deallocate(summ)
             call stop_timer(tmr_std_allocation)
          end do
          if (inode == ionode .and. iprint_basis > 3) &
               write (io_lun, fmt='(6x,"Preconditioned !")')
          call my_barrier()
       else
          Psd = search_direction
       end if
       ! Now determine conjugate directions
       gg = dgg
       dgg = dot(length, search_direction, 1, Psd, 1)
       call gsum(dgg)
       if (inode == ionode .and. iprint_basis > 2) &
            write (io_lun, fmt='(6x,"dgg is ",f20.12)') dgg

       if (gg /= zero) then
          gamma = dgg / gg
       else
          gamma = zero
       end if
       if(gamma > gamma_max) then
         if (inode == ionode) &
              write(io_lun,*) ' Warning: CG direction is reset! '
         gamma = zero
       end if
       if (inode == ionode .and. iprint_basis > 2) &
            write (io_lun, fmt='(6x,"Gamma is ",f20.12)') gamma

       ! Construct the actual search direction
       call copy(length, Psd, 1, search_direction, 1)
       call axpy(length, gamma, last_sd, 1, search_direction, 1)
       ! And project perpendicular to electron gradient
       if (.not. flag_diagonalisation) then
          dN_dot_de = dot(length, search_direction, 1, &
                          elec_grad_coeff_array, 1)
          dN_dot_dN = dot(length, elec_grad_coeff_array, 1, &
                          elec_grad_coeff_array, 1)
          call gsum (dN_dot_de)
          call gsum (dN_dot_dN)
          if (inode == ionode .and. iprint_basis > 3) &
               write (io_lun, fmt='(6x,"dN.de, dN.dN ",2f20.12)') &
                     dN_dot_de, dN_dot_dN
          call axpy(length, - (dN_dot_de / dN_dot_dN), &
                    elec_grad_coeff_array, 1, search_direction, 1)
       end if
       sum_0 = dot(length, grad_coeff_array, 1, search_direction, 1)
       call gsum(sum_0)
       if (inode == ionode .and. iprint_basis > 3) &
            write (io_lun, fmt='(6x,"sum_0 is ",f20.12)') sum_0
       call my_barrier()
       ! minimise the energy (approximately) in this direction.

       if (inode == ionode .and. iprint_basis > 2) &
            write (io_lun, fmt='(6x,"Calling minimise")')

       call my_barrier()

       call line_minimise_support(search_direction, length,      &
                                  fixed_potential, vary_mu,      &
                                  n_L_iterations, tolerance,     &
                                  con_tolerance, total_energy_0, &
                                  expected_reduction, last_step)
       if (inode == ionode .AND. iprint_basis > 2) &
            write (io_lun,fmt='(6x,"Returned !")')
       call dump_blip_coeffs(coefficient_array, coeff_array_size, &
                             inode)
       ! Find change in energy for convergence
       diff = total_energy_last - total_energy_0
       dE_blip = diff
       if (abs(diff / total_energy_0) <= energy_tolerance) then
          if (inode == ionode) &
               write (io_lun, 18) total_energy_0 * en_conv, &
                                  en_units(energy_units)
          convergence_flag = .true.
          total_energy_last = total_energy_0
          return
       end if

       ! prepare for next iteration
       if (n_iterations < n_support_iterations) then
          if (.not. flag_diagonalisation) then
             call LNV_matrix_multiply(electrons, energy_in, doK, doM1, &
                                      doM2, dontM3, doM4, dontphi,     &
                                      dontE, mat_M12=matM12, mat_M4=matM4)
          end if
          if(flag_analytic_blip_int) then
             call get_S_matrix(inode, ionode)
          else
             grad_coeff_array = zero        !! TSUYOSHI MIYAZAKI 2007 Mar 29
          end if
          elec_grad_coeff_array = zero   !! TSUYOSHI MIYAZAKI 2007 Mar 29
          call get_blip_gradient(inode, ionode)
          
          ! The density matrix is effectively idempotent during
          ! diagonalisation
          if (flag_diagonalisation) then
             gridfunctions(H_on_atomfns(1))%griddata = zero
          else
             call get_electron_gradient(atomfns,                &
                                        H_on_atomfns(1), inode, &
                                        ionode)
          end if
       end if

       total_energy_last = total_energy_0

    end do ! n_iterations
    
    deallocate(search_direction, last_sd, Psd, STAT=stat)
    if (stat /= 0) &
         call cq_abort("vary_support: Error dealloc mem")
    call reg_dealloc_mem(area_basis, 3*coeff_array_size, type_dbl)

    call stop_timer (tmr_std_basis)

    return

! 1   format(20x,'mu = ',f10.7,'start energy = ',f15.7)
! 2   format(/20x,'Current Total Energy : ',f15.7,' a.u. ')
! 3   format(20x,'Previous Total Energy: ',f15.7,' a.u. ')
! 4   format(20x,'Difference           : ',f15.7,' a.u. ')
! 5   format(20x,'Required difference  : ',f15.7,' a.u. '/)
7   format(/20x,'------------ Support Variation #: ',i5,' ------------',/)
18  format(///20x,'The blip minimisation has converged to a total energy:', &
         //20x,' Total energy = ',f15.7,' ',a2)

  end subroutine vary_support
  !!***


  !!****f* blip_minimisation/line_minimise_support *
  !!
  !!  NAME 
  !!   line_minimise_support
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Performs a line minimisation on the support functions
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !!   atoms, blip_grid_transform_module, blip,
  !!   calc_matrix_elements_module, common, datatypes, DiagModule,
  !!   dimens, GenBlas, GenComms, logicals, matrix_data,
  !!   maxima_module, mult_module, numbers, SelfCon,
  !!   set_bucket_module, S_matrix_module
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   13/01/98
  !!  MODIFICATION HISTORY
  !!   18/05/2001 dave
  !!    ROBODoc header, changed new_SC_potl call
  !!   24/05/2001 dave
  !!    Shortened call to get_blip_gradient
  !!    Shortened subroutine call
  !!   25/05/2001 dave
  !!    Used get_S_matrix from S_matrix_module
  !!   11/06/2001 dave
  !!    Added RCS Id and Log tags and GenComms dependencies
  !!   17/06/2002 dave
  !!    Added flag to only get K if OrderN solution method is used
  !!    (and tweaked headers)
  !!   31/07/2002 dave
  !!    Changed to use data_M12 from matrix_data and not pass to
  !!    subsidiary routines
  !!   13:52, 04/02/2003 drb 
  !!    Further changes related to diagonalisation (where M12 comes
  !!    from etc)
  !!   09:20, 2003/03/24 dave
  !!    Included in blip_minimisation
  !!   08:29, 2003/04/03 dave
  !!    Changed to use blip_gradient for get_blip_gradient and
  !!    get_electron_gradient
  !!   09:14, 2003/04/10 dave
  !!    Completely rewrote in a more transparent way (closely based
  !!    on safemin in move_atoms.module)
  !!   2008/06/10 ast
  !!    Added timers
  !!   2011/12/05 L.Tong
  !!   - Added spin polarisation
  !!   - Removed redundant parameter number_of_bands
  !!   2012/03/23 L.Tong
  !!   - Changed spin implementation
  !!   - removed redundant input parameter real(double) mu
  !!   2012/06/18 L.Tong
  !!   - removed k0, it seemed not used and redundant
  !!   2017/02/23 dave
  !!    - Changing location of diagon flag from DiagModule to global and name to flag_diagonalisation
  !!   2022/12/14 13:03 dave
  !!    Tweak to make initial k3 larger
  !!  SOURCE
  !!
  subroutine line_minimise_support(search_direction, lengthBlip,  &
                                   fixed_potential, vary_mu,      &
                                   n_cg_L_iterations, tolerance,  &
                                   con_tolerance, total_energy_0, &
                                   expected_reduction, last_step)

    use datatypes
    use numbers
    use logicals
    use global_module,       only: iprint_basis, area_basis, nspin,    &
                                   spin_factor, flag_diagonalisation
    use mult_module,         only: LNV_matrix_multiply, matM12, matM4
    use GenBlas,             only: dot, axpy, copy, vdot
    use SelfCon,             only: new_SC_potl
    use dimens,              only: n_my_grid_points
    use S_matrix_module,     only: get_S_matrix
    use GenComms,            only: gsum, my_barrier, cq_abort, inode,  &
                                   ionode
    ! Check on whether we need K found from L or whether we have it
    ! exactly
    !use DiagModule,          only: diagon
    use support_spec_format, only: coefficient_array,                  &
                                   coeff_array_size, grad_coeff_array, &
                                   elec_grad_coeff_array
    use memory_module,       only: reg_alloc_mem, reg_dealloc_mem,     &
                                   type_dbl

    implicit none

    ! Passed variables
    logical      :: vary_mu, fixed_potential, reset_L
    integer      :: lengthBlip, n_cg_L_iterations
    real(double) :: tolerance, con_tolerance
    real(double) :: total_energy_0
    real(double) :: expected_reduction, last_step
    real(double), dimension(:) :: search_direction

    ! Local variables
    integer      :: i, j, iter, stat, spin
    real(double) :: k1, k2, k3, kmin, lambda, tmp
    real(double) :: e0, e1, e2, e3, energy_out
    real(double), save :: dE = zero ! Use this to guess initial step ?
    real(double), dimension(:), allocatable :: data_blip0
    real(double), dimension(nspin) :: electrons, energy_tmp

    allocate(data_blip0(lengthBlip), STAT=stat)
    if (stat /= 0) &
         call cq_abort("line_minimise_support: Error alloc mem: ", lengthBlip)
    call reg_alloc_mem(area_basis, lengthBlip, type_dbl)

    if (inode == ionode .and. iprint_basis > 2) &
         write (io_lun, &
                fmt='(6x,"On entry to blip line_min, dE is ",e12.5)') dE
    
    call copy(lengthBlip, coefficient_array, 1, data_blip0, 1)
    ! We're assuming that we've ALREADY gone to a self-consistent
    ! ground state before arriving here
    iter = 1
    k1 = zero
    e1 = total_energy_0
    k2 = zero
    e2 = total_energy_0
    e3 = e2
    lambda = two
    ! Loop to find a bracketing triplet
    if (inode == ionode .and. iprint_basis > 0) then
       write (io_lun, fmt='(6x,"Seeking bracketing triplet of points")')
    else if (inode == ionode .and. iprint_basis == 0) then
       write (io_lun, fmt='(6x,"Bracketing")')
    end if
    do while (e3 <= e2)
       if (k2 == zero) then
          k3 = 0.001_double
          ! DRB 2004/03/03
          tmp = vdot(lengthBlip, search_direction, 1, grad_coeff_array, 1)
          if (abs(dE) < RD_ERR) then
             k3 = 0.2_double ! A little conservative
             dE = tmp * k3
          else
             k3 = half * dE / tmp
          end if
          ! elseif (k2==0.01_double) then
          !   k3 = 0.01_double
       else
          k3 = lambda * k2          
       end if
       ! Change blips: start from blip0
       call copy(lengthBlip, data_blip0, 1, coefficient_array, 1)
       call axpy(lengthBlip, k3, search_direction, 1, coefficient_array, 1)
       ! Find new self-consistent energy 
       ! 1. Get new S matrix (includes blip-to-grid transform)
       call get_S_matrix(inode, ionode)
       ! 2. If we're building K as 3LSL-2LSLSL, we need to make K now
       if (.not. flag_diagonalisation) then
          call LNV_matrix_multiply(electrons, energy_tmp, doK, dontM1,&
                                   dontM2, dontM3, dontM4, dontphi,   &
                                   dontE)
       end if
       reset_L = .true.
       ! 3. Get a new self-consistent potential and Hamiltonian
       ! I've not put a call to get_H_matrix here because it's
       ! currently in new_SC_potl
       call new_SC_potl(.false., con_tolerance, reset_L,             &
                        fixed_potential, vary_mu, n_cg_L_iterations, &
                        tolerance, e3)
       if (inode == ionode .and. iprint_basis > 2) &
            write (io_lun, &
                   fmt='(6x,"In blip_min, iter ",i3," step &
                        &and energy are ",2f12.6)') &
                  iter, k3, e3
       if (inode == ionode .and. iprint_basis > 2) &
            write (io_lun, &
                   fmt='(6x," iter=", i3," k1, k2, k3, &
                        &= ",3f12.6)') &
                  iter, k1, k2, k3
       if (e3 < e2) then ! We're still going down hill
          k1 = k2
          e1 = e2
          k2 = k3
          e2 = e3
          iter = iter + 1
       else if (k2 == zero) then
          dE = half * dE
          e3 = e2
       end if
    end do
    ! Interpolate to find minimum.
    kmin = half * (((k1 * k1 - k3 * k3) * (e1 - e2) - &
                    (k1 * k1 - k2 * k2) * (e1 - e3)) / &
                   ((k1 - k3) * (e1 - e2) - (k1 - k2) * (e1 - e3)))
    if (inode == ionode .and. iprint_basis > 1) &
         write (io_lun, &
                fmt='(6x,"In blip_min, bracketed - min from &
                     &extrap: ",4f12.6)') &
               k1, k2, k3, kmin
    if (inode == ionode .and. iprint_basis > 0) then
       write (io_lun, &
              fmt='(6x,"In blip_min, bracketed - energies: ",3f15.6)') &
             e1, e2, e3
    else if (inode == ionode .and. iprint_basis == 0) then
       write (io_lun, fmt='(6x,"Interpolating minimum")')
    end if
    ! Change blips: start from blip0
    call copy(lengthBlip, data_blip0, 1, coefficient_array, 1)
    call axpy(lengthBlip, kmin, search_direction, 1, coefficient_array, 1)
    ! Find new self-consistent energy 
    ! 1. Get new S matrix (includes blip-to-grid transform)
    call get_S_matrix(inode, ionode)
    ! 2. If we're building K as 3LSL-2LSLSL, we need to make K now
    if (.not. flag_diagonalisation) then
       call LNV_matrix_multiply(electrons, energy_tmp, doK, dontM1, &
                                dontM2, dontM3, dontM4, dontphi,    &
                                dontE)
    end if
    reset_L = .true.
    ! 3. Get a new self-consistent potential and Hamiltonian
    ! I've not put a call to get_H_matrix here because it's currently
    ! in new_SC_potl
    call new_SC_potl(.false., con_tolerance, reset_L, fixed_potential,&
                     vary_mu, n_cg_L_iterations, tolerance, energy_out)
    ! If the interpolation failed, go back to the previous "minimum"
    if (energy_out > e2) then 
       kmin = k2
       ! Change blips: start from blip0
       call copy(lengthBlip, data_blip0, 1, coefficient_array, 1)
       call axpy(lengthBlip, kmin, search_direction, 1, &
                 coefficient_array, 1)
       ! Find new self-consistent energy 
       ! 1. Get new S matrix (includes blip-to-grid transform)
       call get_S_matrix(inode, ionode)
       ! 2. If we're building K as 3LSL-2LSLSL, we need to make K now
       if (.not. flag_diagonalisation) then
          call LNV_matrix_multiply(electrons, energy_tmp, doK, dontM1,&
                                   dontM2, dontM3, dontM4, dontphi,   &
                                   dontE)
       end if
       reset_L = .true.
       ! 3. Get a new self-consistent potential and Hamiltonian
       ! I've not put a call to get_H_matrix here because it's
       ! currently in new_SC_potl
       call new_SC_potl(.false., con_tolerance, reset_L,             &
                        fixed_potential, vary_mu, n_cg_L_iterations, &
                        tolerance, energy_out)
    end if
    if (inode == ionode .and. iprint_basis > 0) then
       write (io_lun, fmt='(6x,"In blip_min, at exit energy is ",f15.10)') &
             energy_out
    else if (inode == ionode .and. iprint_basis == 0) then
       write (io_lun, fmt='(6x,"Final Energy: ",f15.10)') energy_out
    end if
    dE = total_energy_0 - energy_out
    if (inode == ionode .and. iprint_basis > 2) &
         write (io_lun, &
                fmt='(6x,"On exit from blip line_min, dE is ",f15.10)') &
               dE
    total_energy_0 = energy_out
    
    deallocate(data_blip0, STAT=stat)
    if (stat /= 0) &
         call cq_abort("line_minimise_support: Error dealloc mem")
    call reg_dealloc_mem(area_basis, lengthBlip, type_dbl)

    return
  end subroutine line_minimise_support
  !!***

end module blip_minimisation
