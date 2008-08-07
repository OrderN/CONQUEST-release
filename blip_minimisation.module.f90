! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
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
!!  SOURCE
!!
module blip_minimisation

  use global_module, ONLY: io_lun
  use timer_stdclocks_module, ONLY: start_timer,stop_timer,tmr_std_basis,tmr_std_allocation

  implicit none

  ! RCS tag for object file identification
  character(len=80), save, private :: RCSid = "$Id$"
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
!!    Changed to use blip_gradient for get_blip_gradient and get_electron_gradient
!!   11:30, 12/11/2004 dave 
!!    Changed to use mx_at_prim from maxima
!!   2007/03/29 08:20 dave and tsuyoshi
!!   2007/04/17 09:36 dave
!!    Added n_L_iterations
!!   2008/06/10 ast
!!    Added timers
!!  SOURCE
!!
  subroutine vary_support( n_support_iterations, fixed_potential, vary_mu, n_L_iterations, &
       number_of_bands, L_tolerance, sc_tolerance, energy_tolerance, mu, total_energy_last, &
       expected_reduction)

    use datatypes
    use logicals
    use numbers
    use units
    use mult_module, ONLY: LNV_matrix_multiply, matM12, matM4
    use GenBlas, ONLY: dot, copy, axpy
    use PosTan, ONLY: PulayC, PulayBeta, SCC, SCBeta
    use blip, ONLY: PreCond, NBlipsRegion
    use GenComms, ONLY: my_barrier, gsum, inode, ionode, cq_abort
    use DiagModule, ONLY: diagon
    use blip_gradient, ONLY: get_blip_gradient, get_electron_gradient
    use global_module, ONLY: flag_precondition_blips, iprint_basis, area_basis
    use io_module, ONLY: dump_blip_coeffs
    use functions_on_grid, ONLY: supportfns, H_on_supportfns, gridfunctions, fn_on_grid
    use support_spec_format, ONLY: coefficient_array, coeff_array_size, grad_coeff_array, elec_grad_coeff_array, &
         supports_on_atom
    use primary_module, ONLY: bundle
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    !     Shared variables

    logical :: vary_mu, fixed_potential
    integer :: n_L_iterations
    real(double) ::  number_of_bands, mu

    integer :: n_support_iterations

    real(double) :: expected_reduction

    real(double) :: total_energy_last, energy_tolerance, L_tolerance, sc_tolerance

    logical convergence_flag

    !     Local variables

    real(double) ::  tolerance, con_tolerance

    integer :: length, n_iterations, n_tries, offset
    logical :: notredone, reduced

    real(double) :: dgg, gamma, gg, electrons, lambda,  &
         step, step_0, step_1, sum_0, sum_1, diff,  &
         total_energy_0, total_energy_test, &
         energy_in, en_dot_el, el_dot_el, kinetic_energy, e_dot_e
    real(double), parameter :: gamma_max = 6.0_double !! TM 2007.03.29

    real(double), allocatable, dimension(:) :: search_direction
    real(double), allocatable, dimension(:) :: last_sd
    real(double), allocatable, dimension(:) :: Psd
    real(double) :: temp,den_del,del_del, last_step, dN_dot_de, dN_dot_dN
    real(double), allocatable, dimension(:) :: sum
    integer indexy, return_ok, jj, n_blip, k, nx, ny, nz, i, j
    integer n,k2,j2, spec, stat

    call start_timer(tmr_std_basis)
    if(inode==ionode) write(io_lun,fmt='(/6x,"Starting blip-coefficient variation",/)')
    if(inode==ionode.AND.iprint_basis>0) &
         write(io_lun,fmt='(6x,"Performing at most ",i4," iterations with a tolerance of ",e12.5,/)') &
         n_support_iterations, energy_tolerance
    call start_timer(tmr_std_allocation)
    allocate(search_direction(coeff_array_size), last_sd(coeff_array_size), Psd(coeff_array_size),STAT=stat)
    if(stat/=0) call cq_abort("Error allocating search directions in vary_support: ",coeff_array_size)
    call reg_alloc_mem(area_basis,3*coeff_array_size,type_dbl)
    call stop_timer(tmr_std_allocation)
    ! Set tolerances for self-consistency and L minimisation
    con_tolerance = SCC*expected_reduction**SCBeta
    tolerance = PulayC*(0.1_double*expected_reduction)**PulayBeta
    if (con_tolerance<sc_tolerance) con_tolerance = sc_tolerance
    if (con_tolerance<10.0_double*tolerance) tolerance = 0.1_double*con_tolerance
    con_tolerance = sc_tolerance
    tolerance = L_tolerance

    if(inode==ionode.AND.iprint_basis>1) write(io_lun,fmt='(6x,"L, SC tolerances: ",2f20.12)') con_tolerance, tolerance

    length = coeff_array_size

    search_direction = zero
    Psd = zero
    grad_coeff_array = zero
    elec_grad_coeff_array = zero

    !     evaluate the action of H on the support functions, calculate
    !     the corresponding matrix elements in the 'support representation'
    !     and calculate the start-point energy

    total_energy_0 = total_energy_last
    if(total_energy_last.eq.zero) total_energy_0 = expected_reduction
    call my_barrier()
    total_energy_last = total_energy_0

    !     now obtain the gradient of the energy with respect to support functions
    ! For diagonalisation, this is already stored
    if(inode==ionode.AND.iprint_basis>2) write(io_lun,fmt='(6x,"Finding gradient")')
    if(.NOT.diagon) call LNV_matrix_multiply(electrons, energy_in, &
         doK, doM1, doM2, dontM3, doM4, dontphi, dontE,matM12,0,matM4,0)
    call get_blip_gradient(inode, ionode)
    if (inode==ionode.AND.iprint_basis>2) write(io_lun,fmt='(6x,"Got blip gradient")')
    call my_barrier()

    ! The density matrix is effectively idempotent during diagonalisation
    if(diagon) then
       gridfunctions(H_on_supportfns)%griddata = zero
    else
       call get_electron_gradient(supportfns, H_on_supportfns, inode, ionode)
       call my_barrier()
       if (inode==ionode.AND.iprint_basis>2) write (*,fmt='(6x,"Got electron gradient")')
    end if

    ! Now we have a basic gradient, so loop
    dgg = zero
    last_step = 1.0e10_double

    !     now loop over search directions

    do n_iterations = 1, n_support_iterations

       if (inode==ionode) then
          if(iprint_basis>0) then
             write(io_lun,7) n_iterations
          else 
             write(io_lun,fmt='(6x,"Support variation: ",i5)') n_iterations
          end if
       end if
       ! We need the last search direction for CG manipulations
       call copy( length, search_direction, 1, last_sd, 1 )

       ! The basis for searching is data_dblip
       call copy( length, grad_coeff_array, 1, search_direction, 1 )
       ! Now project dblip tangential to the constant Ne hyperplane
       if(.NOT.diagon) then
          dN_dot_de = dot( length, grad_coeff_array, 1, elec_grad_coeff_array, 1 )
          dN_dot_dN = dot( length, elec_grad_coeff_array, 1, elec_grad_coeff_array, 1 )
          call gsum(dN_dot_de)
          call gsum(dN_dot_dN)
          if(inode==ionode.AND.iprint_basis>2) write(io_lun,fmt='(6x,"dN.de, dN.dN ",2f20.12)') dN_dot_de, dN_dot_dN
          call axpy( length, -(dN_dot_de/dN_dot_dN), elec_grad_coeff_array, 1, &
               search_direction, 1 )      
       end if
       ! Precondition if we need to
       if(flag_precondition_blips) then
          offset = 0
          if(inode==ionode.AND.iprint_basis>3) write(io_lun,fmt='(6x,"Preconditioning !")')
          call my_barrier()
          do i=1,bundle%n_prim
             spec = bundle%species(i)
             call start_timer(tmr_std_allocation)
             allocate(sum(supports_on_atom(i)%nsuppfuncs))
             call stop_timer(tmr_std_allocation)
             do j=1,NBlipsRegion(spec) ! In future, base on species
                sum=zero
                do k=1,NBlipsRegion(spec) ! Again, base on species
                   do n=1,supports_on_atom(i)%nsuppfuncs
                      sum(n) = sum(n) + PreCond(spec)%coeffs(j,k)*search_direction(offset + (n-1)*NBlipsRegion(spec) + k)
                   enddo
                enddo
                do n=1,supports_on_atom(i)%nsuppfuncs
                   Psd(offset + (n-1)*NBlipsRegion(spec) + j) = sum(n)
                enddo
             enddo
             offset = offset + supports_on_atom(i)%nsuppfuncs*NBlipsRegion(spec)
             call start_timer(tmr_std_allocation)
             deallocate(sum)
             call stop_timer(tmr_std_allocation)
          enddo
          if(inode==ionode.AND.iprint_basis>3) write(io_lun,fmt='(6x,"Preconditioned !")')
          call my_barrier()
       else
          Psd=search_direction
       end if
       ! Now determine conjugate directions
       gg = dgg
       dgg = dot( length,search_direction, 1,Psd, 1)
       call gsum(dgg)
       if(inode==ionode.AND.iprint_basis>2) write(io_lun,fmt='(6x,"dgg is ",f20.12)') dgg

       if (gg.ne.zero) then
          gamma = dgg / gg
       else
          gamma = zero
       end if
       if(gamma > gamma_max) then
         if(inode==ionode) write(io_lun,*) ' Warning: CG direction is reset! '
         gamma=zero
       endif
       if(inode==ionode.AND.iprint_basis>2) write(io_lun,fmt='(6x,"Gamma is ",f20.12)') gamma

       ! Construct the actual search direction
       call copy( length, Psd, 1, search_direction, 1)
       call axpy( length, gamma, last_sd, 1, search_direction, 1)
       ! And project perpendicular to electron gradient
       if(.NOT.diagon) then
          dN_dot_de = dot( length, search_direction, 1, elec_grad_coeff_array, 1 )
          dN_dot_dN = dot( length, elec_grad_coeff_array, 1, elec_grad_coeff_array, 1 )
          call gsum(dN_dot_de)
          call gsum(dN_dot_dN)
          if(inode==ionode.AND.iprint_basis>3) write(io_lun,fmt='(6x,"dN.de, dN.dN ",2f20.12)') dN_dot_de, dN_dot_dN
          call axpy( length, -(dN_dot_de/dN_dot_dN), elec_grad_coeff_array, 1, &
               search_direction, 1 )
       end if
       sum_0 = dot( length, grad_coeff_array, 1, search_direction, 1 )
       call gsum(sum_0)
       if(inode==ionode.AND.iprint_basis>3) write(io_lun,fmt='(6x,"sum_0 is ",f20.12)') sum_0
       call my_barrier()
       ! minimise the energy (approximately) in this direction.

       if(inode==ionode.AND.iprint_basis>2) write(io_lun,fmt='(6x,"Calling minimise")')
       call my_barrier()
       call line_minimise_support( search_direction, length, fixed_potential, vary_mu, n_L_iterations, &
            number_of_bands, tolerance, con_tolerance, mu, total_energy_0, expected_reduction, last_step)
       if(inode==ionode.AND.iprint_basis>2) write(io_lun,'(6x,"Returned !")')
       call dump_blip_coeffs(coefficient_array,coeff_array_size,inode)
       ! Find change in energy for convergence
       diff = total_energy_last - total_energy_0
       if (abs(diff/total_energy_0) .le. energy_tolerance) then
          if (inode==ionode) write(io_lun,18) total_energy_0*en_conv,en_units(energy_units)
          convergence_flag = .true.
          total_energy_last = total_energy_0
          return
       end if

       ! return_ok tells us how the line minimisation went.

       ! prepare for next iteration

       if(.NOT.diagon) call LNV_matrix_multiply(electrons, energy_in, &
            doK, doM1, doM2, dontM3, doM4, dontphi, dontE,matM12,0,matM4,0)
       grad_coeff_array = zero   !! TSUYOSHI MIYAZAKI 2007 Mar 29
       elec_grad_coeff_array = zero   !! TSUYOSHI MIYAZAKI 2007 Mar 29
       call get_blip_gradient(inode, ionode)

       ! The density matrix is effectively idempotent during diagonalisation
       if(diagon) then
          gridfunctions(H_on_supportfns)%griddata = zero
       else
          call get_electron_gradient(supportfns, H_on_supportfns, inode, ionode)
       end if

       total_energy_last = total_energy_0

    end do
    call start_timer(tmr_std_allocation)
    deallocate(search_direction, last_sd, Psd,STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating search directions in vary_support: ",coeff_array_size)
    call reg_dealloc_mem(area_basis,3*coeff_array_size,type_dbl)
    call stop_timer(tmr_std_allocation)
    call stop_timer(tmr_std_basis)

1   format(20x,'mu = ',f10.7,'start energy = ',f15.7)
2   format(/20x,'Current Total Energy : ',f15.7,' a.u. ')
3   format(20x,'Previous Total Energy: ',f15.7,' a.u. ')
4   format(20x,'Difference           : ',f15.7,' a.u. ')
5   format(20x,'Required difference  : ',f15.7,' a.u. '/)
7   format(/20x,'------------ Support Variation #: ',i5,' ------------',/)
18  format(///20x,'The blip minimisation has converged to a total energy:', &
         //20x,' Total energy = ',f15.7,' ',a2)

    return
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
!!   atoms, blip_grid_transform_module, blip, calc_matrix_elements_module, common,
!!   datatypes, DiagModule, dimens, GenBlas, GenComms, logicals, matrix_data, maxima_module,
!!   mult_module, numbers, SelfCon, set_bucket_module, S_matrix_module
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
!!    Added flag to only get K if OrderN solution method is used (and tweaked headers)
!!   31/07/2002 dave
!!    Changed to use data_M12 from matrix_data and not pass to subsidiary routines
!!   13:52, 04/02/2003 drb 
!!    Further changes related to diagonalisation (where M12 comes from etc)
!!   09:20, 2003/03/24 dave
!!    Included in blip_minimisation
!!   08:29, 2003/04/03 dave
!!    Changed to use blip_gradient for get_blip_gradient and get_electron_gradient
!!   09:14, 2003/04/10 dave
!!    Completely rewrote in a more transparent way (closely based on safemin in move_atoms.module)
!!   2008/06/10 ast
!!    Added timers
!!  SOURCE
!!
  subroutine line_minimise_support( search_direction, lengthBlip, fixed_potential, vary_mu, n_cg_L_iterations, &
       number_of_bands, tolerance, con_tolerance, mu, total_energy_0, expected_reduction, last_step )

    use datatypes
    use numbers
    use logicals
    use global_module, ONLY: iprint_basis, area_basis
    use mult_module, ONLY: LNV_matrix_multiply, matM12, matM4
    use GenBlas, ONLY: dot, axpy, copy, vdot
    use SelfCon, ONLY: new_SC_potl
    use dimens, ONLY: n_my_grid_points
    use S_matrix_module, ONLY: get_S_matrix
    use GenComms, ONLY: gsum, my_barrier, cq_abort, inode, ionode
    ! Check on whether we need K found from L or whether we have it exactly
    use DiagModule, ONLY: diagon
    use support_spec_format, ONLY: coefficient_array, coeff_array_size, grad_coeff_array, elec_grad_coeff_array
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    ! Passed variables

    logical :: vary_mu, fixed_potential, reset_L

    integer :: lengthBlip
    integer :: n_cg_L_iterations

    real(double) :: number_of_bands, tolerance, con_tolerance, mu
    real(double) :: total_energy_0
    real(double) :: expected_reduction, last_step
    real(double), dimension(lengthBlip) :: search_direction

    ! Local variables

    real(double), save :: dE = zero ! Use this to guess initial step ?
    real(double), allocatable, dimension(:) :: data_blip0
    real(double) :: k0, k1, k2, k3, kmin, lambda
    real(double) :: e0, e1, e2, e3, electrons, tmp, energy_out
    integer :: i,j, iter, stat

    if(inode==ionode.AND.iprint_basis>2) write(io_lun,fmt='(6x,"On entry to blip line_min, dE is ",e12.5)') dE
    call start_timer(tmr_std_allocation)
    allocate(data_blip0(lengthBlip), STAT=stat)
    if(stat/=0) call cq_abort("LinMinBlip: Error allocating blip storage: ",lengthBlip)
    call reg_alloc_mem(area_basis,lengthBlip,type_dbl)
    call stop_timer(tmr_std_allocation)
    call copy( lengthBlip, coefficient_array, 1, data_blip0, 1)
    ! We're assuming that we've ALREADY gone to a self-consistent ground state before arriving here
    iter = 1
    k1 = zero
    e1 = total_energy_0
    k2 = zero
    e2 = total_energy_0
    e3 = e2
    lambda = 2.0_double
    ! Loop to find a bracketing triplet
    if(inode==ionode.AND.iprint_basis>0) then
       write(io_lun,fmt='(6x,"Seeking bracketing triplet of points")')
    else if(inode==ionode.AND.iprint_basis==0) then
       write(io_lun,fmt='(6x,"Bracketing")')
    end if
    do while(e3<=e2)
       if (k2==zero) then
          k3 = 0.001_double
          ! DRB 2004/03/03
          tmp = vdot(lengthBlip,search_direction,1,grad_coeff_array,1)
          if(abs(dE)<very_small) then
             k3 = 0.008_double
             dE = tmp*k3
          else
             k3 = 0.5_double*dE/tmp
          end if
!       elseif (k2==0.01_double) then
!          k3 = 0.01_double
       else
          k3 = lambda*k2          
       endif
       ! Change blips: start from blip0
       call copy( lengthBlip, data_blip0, 1, coefficient_array, 1)
       call axpy( lengthBlip, k3, search_direction, 1, coefficient_array, 1 )
       ! Find new self-consistent energy 
       ! 1. Get new S matrix (includes blip-to-grid transform)
       call get_S_matrix(inode, ionode)
       ! 2. If we're building K as 3LSL-2LSLSL, we need to make K now
       if(.NOT.diagon) call LNV_matrix_multiply(electrons, tmp, &
            doK, dontM1, dontM2, dontM3, dontM4, dontphi, dontE,0,0,0,0)
       reset_L = .true.
       ! 3. Get a new self-consistent potential and Hamiltonian
       ! I've not put a call to get_H_matrix here because it's currently in new_SC_potl
       call new_SC_potl( .false., con_tolerance, reset_L, fixed_potential, vary_mu, n_cg_L_iterations, &
            number_of_bands, tolerance, mu, e3)
       if(inode==ionode.AND.iprint_basis>2) &
            write(io_lun,fmt='(6x,"In blip_min, iter ",i3," step and energy are ",2f12.6)') iter,k3,e3
       if(inode==ionode.AND.iprint_basis>2) write(io_lun,fmt='(6x," iter=", i3," k0,k1,k2,k3,kmin = ",5f12.6)') &
            iter,k0,k1,k2,k3,kmin
       if (e3<e2) then ! We're still going down hill
          k1 = k2
          e1 = e2
          k2 = k3
          e2 = e3
          iter=iter+1
       else if(k2==zero) then
          dE = 0.5_double*dE
          e3 = e2
       endif
    end do
    ! Interpolate to find minimum.
    kmin = 0.5_double*(((k1*k1-k3*k3)*(e1-e2)-(k1*k1-k2*k2)*(e1-e3))/((k1-k3)*(e1-e2)-(k1-k2)*(e1-e3)))
    if(inode==ionode.AND.iprint_basis>1) &
         write(io_lun,fmt='(6x,"In blip_min, bracketed - min from extrap: ",4f12.6)') k1,k2,k3,kmin
    if(inode==ionode.AND.iprint_basis>0) then
       write(io_lun,fmt='(6x,"In blip_min, bracketed - energies: ",3f15.6)') e1,e2,e3
    else if(inode==ionode.AND.iprint_basis==0) then
       write(io_lun,fmt='(6x,"Interpolating minimum")')
    end if
    ! Change blips: start from blip0
    call copy( lengthBlip, data_blip0, 1, coefficient_array, 1)
    call axpy( lengthBlip, kmin, search_direction, 1, coefficient_array, 1 )
    ! Find new self-consistent energy 
    ! 1. Get new S matrix (includes blip-to-grid transform)
    call get_S_matrix(inode, ionode)
    ! 2. If we're building K as 3LSL-2LSLSL, we need to make K now
    if(.NOT.diagon) call LNV_matrix_multiply(electrons, tmp, &
         doK, dontM1, dontM2, dontM3, dontM4, dontphi, dontE,0,0,0,0)
    reset_L = .true.
    ! 3. Get a new self-consistent potential and Hamiltonian
    ! I've not put a call to get_H_matrix here because it's currently in new_SC_potl
    call new_SC_potl( .false., con_tolerance, reset_L, fixed_potential, vary_mu, n_cg_L_iterations, &
         number_of_bands, tolerance, mu, energy_out)
    ! If the interpolation failed, go back to the previous "minimum"
    if (energy_out>e2) then 
       kmin = k2
       ! Change blips: start from blip0
       call copy( lengthBlip, data_blip0, 1, coefficient_array, 1)
       call axpy( lengthBlip, kmin, search_direction, 1, coefficient_array, 1 )
       ! Find new self-consistent energy 
       ! 1. Get new S matrix (includes blip-to-grid transform)
       call get_S_matrix(inode, ionode)
       ! 2. If we're building K as 3LSL-2LSLSL, we need to make K now
       if(.NOT.diagon) call LNV_matrix_multiply(electrons, tmp, &
            doK, dontM1, dontM2, dontM3, dontM4, dontphi, dontE,0,0,0,0)
       reset_L = .true.
       ! 3. Get a new self-consistent potential and Hamiltonian
       ! I've not put a call to get_H_matrix here because it's currently in new_SC_potl
       call new_SC_potl( .false., con_tolerance, reset_L, fixed_potential, vary_mu, n_cg_L_iterations, &
            number_of_bands, tolerance, mu, energy_out)
    end if
    if(inode==ionode.AND.iprint_basis>0) then
       write(io_lun,fmt='(6x,"In blip_min, at exit energy is ",f15.10)') energy_out
    else if(inode==ionode.AND.iprint_basis==0) then
       write(io_lun,fmt='(6x,"Final Energy: ",f15.10)') energy_out
    endif
    dE = total_energy_0 - energy_out
    if(inode==ionode.AND.iprint_basis>2) write(io_lun,fmt='(6x,"On exit from blip line_min, dE is ",f15.10)') dE
    total_energy_0 = energy_out
    call start_timer(tmr_std_allocation)
    deallocate(data_blip0, STAT=stat)
    if(stat/=0) call cq_abort("LinMinBlip: Error allocating blip storage: ",lengthBlip)
    call reg_dealloc_mem(area_basis,lengthBlip,type_dbl)
    call stop_timer(tmr_std_allocation)
    return
  end subroutine line_minimise_support
!!***

end module blip_minimisation
