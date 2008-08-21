! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module pao_minimisation
! ------------------------------------------------------------------------------
! Code area 6: energy minimisation
! ------------------------------------------------------------------------------

!!****h* Conquest/pao_minimisation *
!!  NAME
!!   pao_minimisation
!!  PURPOSE
!!   Hold the different routines associated with pao minimisation
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   09:17, 2003/03/24 dave
!!  MODIFICATION HISTORY
!!   2005/07/11 10:22 dave
!!    Changed build_PAO_coeff_grad so that nsf and npao are passed in with 
!!    sensible names and tidied module use
!!   10:09, 13/02/2006 drb 
!!    Removed all explicit references to data_ variables and rewrote in terms of new 
!!    matrix routines
!!   2008/02/06 08:31 dave
!!    Changed for output to file not stdout
!!   2008/05/25 ast
!!    Added timers
!!  SOURCE
!!
module pao_minimisation

  use datatypes
  use global_module, ONLY: io_lun
  use timer_stdclocks_module, ONLY: start_timer,stop_timer,tmr_std_allocation,tmr_std_matrices

  implicit none

  integer, parameter :: mx_pulay = 5
  integer, parameter :: GdS = 1
  integer, parameter :: KdH = 2
  integer, parameter :: full = 3
  !real(double), save  :: PAOprecond(npao,npao,mx_at_prim*mx_tnab)
  !real(double), save  :: PAOprecond2(npao)

  ! RCS tag for object file identification
  character(len=80), save, private :: RCSid = "$Id$"
!!***

contains

!!****f* pao_minimisation/vary_pao *
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
!!   pao, common, datatypes, DiagModule, GenBlas, GenComms, logicals, 
!!   matrix_data, maxima_module, mult_module, numbers, PosTan
!!  AUTHOR
!!   D.R.Bowler/E.H.Hernandez
!!  CREATION DATE
!!   03/04/95
!!  MODIFICATION HISTORY
!!   updated to make use of pao functions.
!!   updated to include electron gradient
!!   25/1/96 CMG/EHE line minimisation moved out
!!   19/5/97 DRB added HeadGordon
!!   23/05/2001 dave
!!    Shortened calls to get_pao_gradient and get_electron_gradient
!!    and added ROBODoc header
!!   23/05/2001 dave
!!    Shortened call to line_minimise_support
!!   08/06/2001 dave
!!    Added RCS Id and Log tags and changed to use GenComms
!!   17/06/2002 dave
!!    Tidied headers, added check for solution method
!!   09:18, 2003/03/24 dave
!!    Included in pao_minimisation
!!   08:29, 2003/04/03 dave
!!    Changed to use pao_gradient for get_pao_gradient and get_electron_gradient
!!   2007/04/26 12:08 dave
!!    Changed TestPAOGrads to TestBasisGrads (to allow both blip and PAO testing with same flag)
!!   2008/05/25 ast
!!    Added timers
!!  SOURCE
!!
  subroutine vary_pao( n_support_iterations, fixed_potential, vary_mu, n_cg_L_iterations, &
       number_of_bands, L_tolerance, sc_tolerance, energy_tolerance, mu, total_energy_last, expected_reduction)

    use datatypes
    use logicals
    use numbers
    use GenBlas
    use PosTan, ONLY: PulayC, PulayBeta, SCC, SCBeta
    use GenComms, ONLY: my_barrier, gsum, inode, ionode, cq_abort
    use DiagModule, ONLY: diagon
    use primary_module, ONLY: bundle
    use global_module, ONLY: flag_vary_basis, iprint_minE, ni_in_cell, flag_self_consistent, id_glob, numprocs
    use group_module, ONLY: parts
    use H_matrix_module, ONLY: get_H_matrix
    use S_matrix_module, ONLY: get_S_matrix
    use make_rad_tables, ONLY: writeout_support_functions
    use support_spec_format, ONLY: supports_on_atom, coeff_array_size, grad_coeff_array, elec_grad_coeff_array,&
         support_gradient, support_elec_gradient, flag_paos_atoms_in_cell, &
         TestBasisGrads, TestTot, TestBoth, TestS, TestH
    use DMMin, ONLY: FindMinDM
    use energy, ONLY: get_energy, kinetic_energy, nl_energy, band_energy
    use PAO_grid_transform_module, ONLY: PAO_to_grid
    use functions_on_grid, ONLY: supportfns
    use species_module, ONLY: nsf_species, npao_species
    use density_module, ONLY: density
    use maxima_module, ONLY: maxngrid

    implicit none

    !     Shared variables
    logical vary_mu, fixed_potential
    integer n_cg_L_iterations
    real(double) ::  number_of_bands, mu

    integer :: n_support_iterations
    real(double) :: expected_reduction

    real(double) :: total_energy_last, energy_tolerance, L_tolerance, sc_tolerance
    logical convergence_flag

    !     Local variables
    real(double) ::  tolerance, con_tolerance, tmp

    integer :: length, n_iterations, n_tries
    logical :: notredone, reduced, orig_SC, reset_L

    real(double) :: dgg, gamma, gg, electrons, lambda,  &
         step, step_0, step_1, sum_0, sum_1, diff,  &
         total_energy_0, total_energy_test, &
         energy_in, en_dot_el, el_dot_el, e_dot_e

    real(double), dimension(:), allocatable :: search_direction, last_sd, Psd
    real(double), dimension(:), allocatable :: grad_copy, grad_copy_dH, grad_copy_dS
    real(double) :: temp,den_del,del_del, last_step, dN_dot_de, dN_dot_dN
    real(double) :: sum, tmpgrad, E2, E1, H1, H2, H1a, H2a, BE2, BE1, g1, g2

    integer indexy, return_ok, jj, n_blip, k, nx, ny, nz, i, j
    integer n,k2,j2, part, memb, nsf1, npao1, point, iprim, nab, nsf2, proc, ind_part, atom, which_atom, local_atom

    logical :: my_atom

    reset_L = .true.
    call start_timer(tmr_std_allocation)
    allocate(search_direction(coeff_array_size),last_sd(coeff_array_size), Psd(coeff_array_size))
    if(TestBasisGrads) then
       allocate(grad_copy(coeff_array_size),grad_copy_dH(coeff_array_size),grad_copy_dS(coeff_array_size))
    end if
    call stop_timer(tmr_std_allocation)
    ! Set tolerances for self-consistency and L minimisation
    con_tolerance = zero ! SCC*expected_reduction**SCBeta
    tolerance = zero ! PulayC*(0.1_double*expected_reduction)**PulayBeta
    if (con_tolerance<sc_tolerance) con_tolerance = sc_tolerance
    if (con_tolerance<10.0_double*tolerance) tolerance = 0.1_double*con_tolerance
    con_tolerance = sc_tolerance

    if(inode==ionode) write(io_lun,*) 'Tolerances: ',con_tolerance, tolerance
    if (inode.eq.ionode) write(io_lun,*) INODE,' entering vary_pao'

    length = coeff_array_size

    search_direction = 0.0_double
    Psd = 0.0_double
    last_sd = 0.0_double ! TO

    total_energy_0 = total_energy_last
    if(total_energy_last.eq.0.0_double) total_energy_0 = expected_reduction
    total_energy_last = total_energy_0

    ! We need to assemble the gradient
    grad_coeff_array = zero
    elec_grad_coeff_array = zero
    ! call get _H_matrix before calling build_PAO_coeff ! TM
    call get_H_matrix(.true., fixed_potential, electrons, density, maxngrid)
    call build_PAO_coeff_grad(full)
    if(TestBasisGrads) then
       grad_copy = grad_coeff_array
       do i = 1,size(grad_coeff_array)
          write(20,*) grad_coeff_array(i)
       end do
       ! Test PAO gradients
       ! Preserve unperturbed energy and gradient
       E1 = band_energy
       grad_coeff_array = zero
       elec_grad_coeff_array = zero
       call build_PAO_coeff_grad(GdS)
       do i = 1,size(grad_coeff_array)
          write(21,*) grad_coeff_array(i)
       end do
       grad_copy_dS = grad_coeff_array
       grad_coeff_array = zero
       elec_grad_coeff_array = zero
       call build_PAO_coeff_grad(KdH)
       do i = 1,size(grad_coeff_array)
          write(22,*) grad_coeff_array(i)
       end do
       grad_copy_dH = grad_coeff_array
       do proc = 1, numprocs
          local_atom = 0
          if(inode==proc) then 
             my_atom = .true.
          else
             my_atom = .false.
          end if
          do part = 1, parts%ng_on_node(proc)
             ind_part = parts%ngnode(parts%inode_beg(proc)+part-1)
             do atom = 1, parts%nm_group(ind_part)
                ! Need global number and local number
                local_atom = local_atom+1 ! Is this really primary atom ?
                if(flag_paos_atoms_in_cell) then
                   which_atom = id_glob(parts%icell_beg(ind_part)+atom-1)
                else
                   which_atom = local_atom
                end if
                if(my_atom) write(io_lun,*) 'global, primary, prim(glob) ',which_atom, local_atom, bundle%ig_prim(local_atom)
                do nsf1 = 1,nsf_species(bundle%species(local_atom))
                   do npao1 = 1,npao_species(bundle%species(local_atom))
                      tmp = 0.0001_double*supports_on_atom(which_atom)%supp_func(nsf1)%coefficients(npao1)
                      if(TestTot) then
                         ! Shift coefficient a little
                         supports_on_atom(which_atom)%supp_func(nsf1)%coefficients(npao1) = &
                              supports_on_atom(which_atom)%supp_func(nsf1)%coefficients(npao1) + tmp
                         call my_barrier()
                         ! Recalculate energy and gradient
                         call get_S_matrix(inode, ionode)
                         call get_H_matrix(.true., fixed_potential, electrons, density, maxngrid)
                         call FindMinDM(n_cg_L_iterations, number_of_bands, vary_mu, &
                              L_tolerance, mu, inode, ionode, .false., .false.)
                         call get_energy(E2)
                         E2 = band_energy
                         grad_coeff_array = zero
                         elec_grad_coeff_array = zero
                         call build_PAO_coeff_grad(full)
                         if(my_atom) then
                            g1 = support_gradient(which_atom)%supp_func(nsf1)%coefficients(npao1)
                            grad_coeff_array = grad_copy
                            g2 = support_gradient(which_atom)%supp_func(nsf1)%coefficients(npao1)
                            write(io_lun,*) 'Tot: Numerical, analytic grad: ', (E2-E1)/tmp, -0.5_double*(g1+g2)
                            write(io_lun,*) 'Tot:Components: ', tmp,E1,E2,g1,g2
                         end if
                         ! Shift coefficient back
                         supports_on_atom(which_atom)%supp_func(nsf1)%coefficients(npao1) = &
                              supports_on_atom(which_atom)%supp_func(nsf1)%coefficients(npao1) - tmp
                         call get_S_matrix(inode, ionode)
                         call get_H_matrix(.true., fixed_potential, electrons, density, maxngrid)
                         call my_barrier()
                      end if
                      if(TestS.or.TestBoth) then
                         ! Shift coefficient a little
                         !tmp = 0.0001_double*supports_on_atom(which_atom)%supp_func(nsf1)%coefficients(npao1)
                         supports_on_atom(which_atom)%supp_func(nsf1)%coefficients(npao1) = &
                              supports_on_atom(which_atom)%supp_func(nsf1)%coefficients(npao1) + tmp
                         call my_barrier()
                         ! Recalculate energy and gradient
                         call get_S_matrix(inode, ionode)
                         call FindMinDM(n_cg_L_iterations, number_of_bands, vary_mu, &
                              L_tolerance, mu, inode, ionode, .false., .false.)
                         call get_energy(E2)
                         E2 = band_energy
                         grad_coeff_array = zero
                         elec_grad_coeff_array = zero
                         call build_PAO_coeff_grad(GdS)
                         if(my_atom) then
                            g1 = support_gradient(which_atom)%supp_func(nsf1)%coefficients(npao1)
                            grad_coeff_array = grad_copy_dS
                            g2 = support_gradient(which_atom)%supp_func(nsf1)%coefficients(npao1)
                            write(io_lun,*) 'GdS: Numerical, analytic grad: ', (E2-E1)/tmp, -0.5_double*(g1+g2)
                            write(io_lun,*) 'GdS:Components: ', tmp,E1,E2,g1,g2
                         end if
                         ! Shift coefficient back
                         supports_on_atom(which_atom)%supp_func(nsf1)%coefficients(npao1) = &
                              supports_on_atom(which_atom)%supp_func(nsf1)%coefficients(npao1) - tmp
                         call get_S_matrix(inode, ionode)
                         call my_barrier()
                      end if
                      ! ** Test H ** !
                      if(TestH.or.TestBoth) then
                         ! Shift coefficient a little
                         !tmp = 0.0001_double*supports_on_atom(which_atom)%supp_func(nsf1)%coefficients(npao1)
                         supports_on_atom(which_atom)%supp_func(nsf1)%coefficients(npao1) = &
                              supports_on_atom(which_atom)%supp_func(nsf1)%coefficients(npao1) + tmp
                         call my_barrier()
                         ! Recalculate energy and gradient
                         call PAO_to_grid(inode-1,supportfns)
                         call get_H_matrix(.true., fixed_potential, electrons, density, maxngrid)
                         call FindMinDM(n_cg_L_iterations, number_of_bands, vary_mu, &
                              L_tolerance, mu, inode, ionode, .false., .false.)
                         call get_energy(E2)
                         E2 = band_energy
                         grad_coeff_array = zero
                         elec_grad_coeff_array = zero
                         call build_PAO_coeff_grad(KdH)
                         if(my_atom) then
                            g1 = support_gradient(which_atom)%supp_func(nsf1)%coefficients(npao1)
                            grad_coeff_array = grad_copy_dH
                            g2 = support_gradient(which_atom)%supp_func(nsf1)%coefficients(npao1)
                            write(io_lun,*) 'KdH: Numerical, analytic grad: ', (E2-E1)/tmp, -0.5_double*(g1+g2)
                            write(io_lun,*) 'KdH:Components: ', tmp,E1,E2,g1,g2
                         end if
                         ! Shift coefficient back
                         supports_on_atom(which_atom)%supp_func(nsf1)%coefficients(npao1) = &
                              supports_on_atom(which_atom)%supp_func(nsf1)%coefficients(npao1) - tmp
                         call PAO_to_grid(inode-1,supportfns)
                         call get_H_matrix(.true., fixed_potential, electrons, density, maxngrid)
                      end if
                   enddo
                enddo
             enddo
          enddo
       enddo
       call start_timer(tmr_std_allocation)
       deallocate(grad_copy,grad_copy_dH,grad_copy_dS)
       call stop_timer(tmr_std_allocation)
    end if ! TestBasisGrads
    ! What about preconditioning ?       
    call my_barrier()
    ! Now we have a basic gradient, so loop
    dgg = zero
    last_step = 1.0D10

    !     now loop over search directions
    do n_iterations = 1, n_support_iterations
       if (inode .eq. ionode) write(io_lun,7) n_iterations
       ! We need the last search direction for CG manipulations
       call copy( length, search_direction, 1, last_sd, 1 )
       ! The basis for searching is gradient
       call copy( length, grad_coeff_array, 1, search_direction, 1 )
       ! Now project gradient tangential to the constant Ne hyperplane
       if(.NOT.diagon) then
          dN_dot_de = dot( length, grad_coeff_array, 1, elec_grad_coeff_array, 1 )
          dN_dot_dN = dot( length, elec_grad_coeff_array, 1, elec_grad_coeff_array, 1 )
          if(.NOT.flag_paos_atoms_in_cell) then
             call gsum(dN_dot_de)
             call gsum(dN_dot_dN)
          end if
          if(inode.eq.ionode) write(io_lun,*) 'dN.de, dN.dN ',dN_dot_de, dN_dot_dN
          call axpy( length, -(dN_dot_de/dN_dot_dN), elec_grad_coeff_array, 1, &
               search_direction, 1 )      
       end if
       ! *THINK* Do we need/want to precondition ?
       Psd=search_direction
       ! Now determine conjugate directions
       gg = dgg
       dgg = dot( length,search_direction, 1,Psd, 1)
       call gsum(dgg)
       if(inode.eq.ionode) write(io_lun,*) 'dgg is ',dgg
       if (gg.ne.zero) then
          gamma = dgg / gg
       else
          gamma = zero
       end if
       !gamma = zero
       !if(mod(n_iterations,5)==0) gamma = zero
       if(inode.eq.ionode) write(io_lun,*) 'Gamma is ',gamma

       ! Construct the actual search direction
       call copy( length, Psd, 1, search_direction, 1)
       call axpy( length, gamma, last_sd, 1, search_direction, 1)
       ! And project perpendicular to electron gradient
       if(.NOT.diagon) then
          dN_dot_de = dot( length, search_direction, 1, elec_grad_coeff_array, 1 )
          dN_dot_dN = dot( length, elec_grad_coeff_array, 1, elec_grad_coeff_array, 1 )
          if(.NOT.flag_paos_atoms_in_cell) then
             call gsum(dN_dot_de)
             call gsum(dN_dot_dN)
          end if
          if(inode.eq.ionode) write(io_lun,*) 'dN.de, dN.dN ',dN_dot_de, dN_dot_dN
          call axpy( length, -(dN_dot_de/dN_dot_dN), elec_grad_coeff_array, 1, &
               search_direction, 1 )
       end if
       ! Check this !
       sum_0 = dot( length, grad_coeff_array, 1, search_direction, 1 )
       if(.NOT.flag_paos_atoms_in_cell) call gsum(sum_0)
       if(inode.eq.ionode) write(io_lun,*) 'sum_0 is ',sum_0
       call my_barrier()

       ! minimise the energy (approximately) in this direction.
       if(inode.eq.ionode) write(io_lun,*) 'Minimise'
       call my_barrier()
       tmp = dot(length,search_direction,1,grad_coeff_array,1)
       if(.NOT.flag_paos_atoms_in_cell) call gsum(tmp)
       ! Temporarily turn off basis variationso that we don't do unnecessary calculations
       flag_vary_basis = .false.
       !orig_SC = flag_self_consistent
       !flag_self_consistent = .false.
       call line_minimise_pao( search_direction, fixed_potential, vary_mu, n_cg_L_iterations, &
            number_of_bands, tolerance, con_tolerance, mu, total_energy_0, expected_reduction, last_step, tmp)
       if(inode.eq.ionode) write(io_lun,*) 'Returned !'
       do i=1,ni_in_cell
          do nsf1 = 1,supports_on_atom(i)%nsuppfuncs
             sum = zero
             do npao1 = 1,supports_on_atom(i)%supp_func(nsf1)%ncoeffs
                sum = sum + supports_on_atom(i)%supp_func(nsf1)%coefficients(npao1)* &
                     supports_on_atom(i)%supp_func(nsf1)%coefficients(npao1)
             end do
             sum = sqrt(sum)
             supports_on_atom(i)%supp_func(nsf1)%coefficients = supports_on_atom(i)%supp_func(nsf1)%coefficients/sum
          end do
       end do
       call writeout_support_functions(inode,ionode)
       ! Find change in energy for convergence
       diff = total_energy_last - total_energy_0
       if (abs(diff/total_energy_0) .le. energy_tolerance) then
          if (inode .eq. ionode) write(io_lun,18) total_energy_0
          convergence_flag = .true.
          total_energy_last = total_energy_0
          return
       end if
       ! prepare for next iteration
       flag_vary_basis = .true.
       ! Find new self-consistent energy 
       ! 1. Generate data_dS
       call get_S_matrix(inode, ionode)
       ! 3. Generate data_dH
       call get_H_matrix(.true., fixed_potential, electrons, density, maxngrid)
       call FindMinDM(n_cg_L_iterations, number_of_bands, vary_mu, &
            L_tolerance, mu, inode, ionode, .false., .false.)
       write(io_lun,*) 'HERE'
       call get_energy(total_energy_test)
       ! We need to assemble the gradient
       grad_coeff_array = zero
       elec_grad_coeff_array = zero
       call build_PAO_coeff_grad(full)
       !if(inode==1) then
       !   gradient(:,:,2:mx_at_prim) = zero
       !else
       !   gradient = zero
       !end if
       sum = dot(length,grad_coeff_array,1,grad_coeff_array,1)
       if(.NOT.flag_paos_atoms_in_cell) call gsum(sum)
       write(io_lun,*) 'Dot prod of gradient: ',sum
       total_energy_last = total_energy_0
    end do

1   format(20x,'mu = ',f10.7,'start energy = ',f15.7)
2   format(/20x,'Current Total Energy : ',f15.7,' a.u. ')
3   format(20x,'Previous Total Energy: ',f15.7,' a.u. ')
4   format(20x,'Difference           : ',f15.7,' a.u. ')
5   format(20x,'Required difference  : ',f15.7,' a.u. '/)
7   format(/20x,'------------ PAO Variation #: ',i5,' ------------',/)
18  format(///20x,'The minimisation has converged to a total energy:', &
         //20x,' Total energy = ',f15.7)

    return
  end subroutine vary_pao
!!***


  subroutine pulay_min_pao( n_support_iterations, fixed_potential, vary_mu, n_cg_L_iterations, &
       number_of_bands, L_tolerance, sc_tolerance, energy_tolerance, mu, total_energy_last, expected_reduction)

    use datatypes
    use logicals
    use numbers
    use mult_module, ONLY: LNV_matrix_multiply, matM12, matM4
    use GenBlas, ONLY: dot, copy
    use PosTan, ONLY: PulayC, PulayBeta, SCC, SCBeta
    use GenComms, ONLY: my_barrier, gsum, inode, ionode
    use DiagModule, ONLY: diagon
    use primary_module, ONLY: bundle
    use global_module, ONLY: flag_vary_basis, iprint_minE, ni_in_cell
    use SelfCon, ONLY: new_SC_potl
    use S_matrix_module, ONLY: get_S_matrix
    use make_rad_tables, ONLY: writeout_support_functions
    use support_spec_format, ONLY: supports_on_atom
    use DMMin, ONLY: FindMinDM
    use energy, ONLY: get_energy, kinetic_energy, nl_energy
    use support_spec_format, ONLY: supports_on_atom, coeff_array_size, grad_coeff_array, elec_grad_coeff_array,&
         support_gradient, support_elec_gradient, coefficient_array, flag_paos_atoms_in_cell
    use Pulay

    implicit none

    !     Shared variables
    logical vary_mu, fixed_potential
    integer n_cg_L_iterations
    real(double) ::  number_of_bands, mu

    integer ::  n_support_iterations
    real(double) :: expected_reduction

    real(double) :: total_energy_last, energy_tolerance, L_tolerance, sc_tolerance
    logical convergence_flag

    !     Local variables
    real(double) ::  tolerance, con_tolerance, tmp

    integer :: length, n_iterations, n_tries, npmod, pul_mx
    logical :: notredone, reduced, reset_L

    real(double) :: dgg, gamma, gg, electrons, lambda,  &
         step, step_0, step_1, sum_0, sum_1, diff,  &
         total_energy_0, total_energy_test, &
         energy_in, en_dot_el, el_dot_el, e_dot_e, deltaE, g0
    real(double) :: Aij(mx_pulay, mx_pulay), alph(mx_pulay)
    real(double) :: Aij1(mx_pulay*mx_pulay)

    real(double), dimension(:), allocatable :: search_direction, last_sd,Psd
    real(double) :: temp,den_del,del_del, last_step, dN_dot_de, dN_dot_dN
    real(double) :: sum, tmpgrad, E2, E1, H1, H2, H1a, H2a
    integer indexy, return_ok, jj, n_blip, k, nx, ny, nz, i, j, ii
    integer n,k2,j2, part, memb, nsf1, npao1, point, iprim, nab, nsf2
    real(double) :: data_gradstore(coeff_array_size,mx_pulay)
    real(double) :: data_paostore(coeff_array_size,mx_pulay)


    call start_timer(tmr_std_allocation)
    allocate(search_direction(coeff_array_size),last_sd(coeff_array_size), Psd(coeff_array_size))
    call stop_timer(tmr_std_allocation)
    ! Set tolerances for self-consistency and L minimisation
    con_tolerance = SCC*expected_reduction**SCBeta
    tolerance = PulayC*(0.1_double*expected_reduction)**PulayBeta
    if (con_tolerance<sc_tolerance) con_tolerance = sc_tolerance
    if (con_tolerance<10.0_double*tolerance) tolerance = 0.1_double*con_tolerance
    con_tolerance = sc_tolerance

    if(inode==ionode) write(io_lun,*) 'Tolerances: ',con_tolerance, tolerance, energy_tolerance
    if (inode.eq.ionode) write(io_lun,*) INODE,' entering vary_pao'

    length = coeff_array_size

    search_direction = zero
    Psd = zero

    total_energy_0 = total_energy_last
    if(total_energy_last.eq.0.0_double) total_energy_0 = expected_reduction
    total_energy_last = total_energy_0

    ! We need to assemble the gradient
    if(.NOT.diagon) call LNV_matrix_multiply(electrons, energy_in, &
         doK, doM1, doM2, dontM3, doM4, dontphi, dontE,matM12,0,matM4,0)
    ! We should have the elements built by H_matrix_module and S_matrix_module
    ! Now we take the sum over j\beta (nsf2 = \beta; neigh = j)
    grad_coeff_array = zero
    elec_grad_coeff_array = zero
    call build_PAO_coeff_grad(full)
    ! What about preconditioning ?       
    call my_barrier()
    ! Now we have a basic gradient, so loop
    g0 = dot(length,grad_coeff_array,1,grad_coeff_array,1)
    if(.NOT.flag_paos_atoms_in_cell) call gsum(g0)
    last_step = 1.0D10
    write(io_lun,*) 'Dot product of initial gradient ',g0
    ! Store gradient
    call copy(length,grad_coeff_array,1,data_gradstore(1:,1),1)
    data_paostore(:,1) = coefficient_array
    diff = zero
    !     now loop over search directions
    do n_iterations = 1, n_support_iterations
       if (inode .eq. ionode) write(io_lun,7) n_iterations
       npmod = mod(n_iterations, mx_pulay)+1
       pul_mx = min(n_iterations+1, mx_pulay)
       step = diff/g0 ! Base step on present gradient and expected dE
       if(step == zero) step = 0.01_double
       if(inode==ionode) write(io_lun,*) 'npmod, pul_mx and step: ',npmod, pul_mx,step
       ! Build PAO coefficients
       if(npmod>1) then
          coefficient_array = coefficient_array + step*data_gradstore(:,npmod-1)
       else
          coefficient_array = coefficient_array + step*data_gradstore(:,pul_mx)
       endif
       write(io_lun,*) 'Normalising'
       ! Normalise
       do ii=1,ni_in_cell
          do nsf1 = 1,supports_on_atom(ii)%nsuppfuncs ! Select alpha
             sum = zero
             do npao1 = 1,supports_on_atom(ii)%supp_func(nsf1)%ncoeffs ! PAOs for i, alpha
                sum = sum + supports_on_atom(ii)%supp_func(nsf1)%coefficients(npao1)* &
                     supports_on_atom(ii)%supp_func(nsf1)%coefficients(npao1)
             end do
             sum = sqrt(sum)
             supports_on_atom(ii)%supp_func(nsf1)%coefficients = supports_on_atom(ii)%supp_func(nsf1)%coefficients/sum
          end do
       end do
       ! Find change in energy for convergence
       ! Get energy and gradient for step
       write(io_lun,*) 'Getting new S, H, E'
       flag_vary_basis = .true.
       ! Find new self-consistent energy 
       ! 1. Generate data_dS
       call get_S_matrix(inode, ionode)
       ! 2. If we're building K as 3LSL-2LSLSL, we need to make K now
       if(.NOT.diagon) call LNV_matrix_multiply(electrons, tmp, &
            doK, dontM1, dontM2, dontM3, dontM4, dontphi, dontE,0,0,0,0)
       reset_L = .true.
       ! 3. Get a new self-consistent potential and Hamiltonian
       ! I've not put a call to get_H_matrix here because it's currently in new_SC_potl
       call new_SC_potl( .false., con_tolerance, reset_L, fixed_potential, vary_mu, n_cg_L_iterations, &
            number_of_bands, tolerance, mu, total_energy_0)
       grad_coeff_array = zero
       elec_grad_coeff_array = zero
       call build_PAO_coeff_grad(full)
       sum = dot(length,grad_coeff_array,1,grad_coeff_array,1)
       if(.NOT.flag_paos_atoms_in_cell) call gsum(sum)
       write(io_lun,*) 'Dot prod of gradient: ',sum
       ! Store PAO and gradient at this step
       call copy(length,grad_coeff_array,1,data_gradstore(1:,npmod),1)
       call copy(length,coefficient_array,1,data_paostore(1:,npmod),1)
       ! Now mix pulay
       Aij = zero
       do ii=1,pul_mx
          do j=1,pul_mx
             gg = dot(length, data_gradstore(1:,j),1, &
                  data_gradstore(1:,ii),1)
             if(.NOT.flag_paos_atoms_in_cell) call gsum(gg)
             Aij(j,ii) = gg
             Aij1(j+(ii-1)*pul_mx) = gg
          enddo
       enddo
       ! Solve to get alphas
       call DoPulay2D(Aij,alph,pul_mx,mx_pulay,inode,ionode)
       write(io_lun,*) 'Alph: ',alph
       ! Make new supports
       coefficient_array = zero
       do ii=1,pul_mx
          coefficient_array = coefficient_array + alph(ii)*data_paostore(:,ii)
       end do
       ! re-evaluate the gradient and energy at new position
       ! Find new self-consistent energy 
       ! 1. Generate data_dS
       call get_S_matrix(inode, ionode)
       ! 2. If we're building K as 3LSL-2LSLSL, we need to make K now
       if(.NOT.diagon) call LNV_matrix_multiply(electrons, tmp, &
            doK, dontM1, dontM2, dontM3, dontM4, dontphi, dontE,0,0,0,0)
       reset_L = .true.
       ! 3. Get a new self-consistent potential and Hamiltonian
       ! I've not put a call to get_H_matrix here because it's currently in new_SC_potl
       call new_SC_potl( .false., con_tolerance, reset_L, fixed_potential, vary_mu, n_cg_L_iterations, &
            number_of_bands, tolerance, mu, total_energy_0)
       ! Normalise
       do ii=1,ni_in_cell
          do nsf1 = 1, supports_on_atom(ii)%nsuppfuncs! Select alpha
             sum = zero
             do npao1 = 1,supports_on_atom(ii)%supp_func(nsf1)%ncoeffs ! PAOs for i, alpha
                sum = sum + supports_on_atom(ii)%supp_func(nsf1)%coefficients(npao1)* &
                     supports_on_atom(ii)%supp_func(nsf1)%coefficients(npao1)
             end do
             sum = sqrt(sum)
             supports_on_atom(ii)%supp_func(nsf1)%coefficients = supports_on_atom(ii)%supp_func(nsf1)%coefficients/sum
          end do
       end do
       grad_coeff_array = zero
       elec_grad_coeff_array = zero
       call build_PAO_coeff_grad(full)
       sum = dot(length,grad_coeff_array,1,grad_coeff_array,1)
       if(.NOT.flag_paos_atoms_in_cell) call gsum(sum)
       write(io_lun,*) 'Dot prod of gradient: ',sum
       ! Replace step with real L
       call copy(length,grad_coeff_array,1,data_gradstore(1:,npmod),1)
       call copy(length,coefficient_array,1,data_paostore(1:,npmod),1)
       call writeout_support_functions(inode,ionode)
       diff = total_energy_last - total_energy_0
       total_energy_last = total_energy_0
       if (abs(diff/total_energy_0) .le. energy_tolerance) then
          if (inode .eq. ionode) write(io_lun,18) total_energy_0
          convergence_flag = .true.
          total_energy_last = total_energy_0
          return
       end if
    end do

1   format(20x,'mu = ',f10.7,'start energy = ',f15.7)
2   format(/20x,'Current Total Energy : ',f15.7,' a.u. ')
3   format(20x,'Previous Total Energy: ',f15.7,' a.u. ')
4   format(20x,'Difference           : ',f15.7,' a.u. ')
5   format(20x,'Required difference  : ',f15.7,' a.u. '/)
7   format(/20x,'------------ PAO Variation #: ',i5,' ------------',/)
18  format(///20x,'The minimisation has converged to a total energy:', &
         //20x,' Total energy = ',f15.7)

    return
  end subroutine pulay_min_pao
!!***


!!****f* pao_minimisation/line_minimise_pao *
!!
!!  NAME 
!!   line_minimise_pao
!!  USAGE
!! 
!!  PURPOSE
!!   Performs a line minimisation on the support functions
!!  INPUTS
!! 
!! 
!!  USES
!!   atoms, pao_grid_transform_module, pao, calc_matrix_elements_module, common,
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
!!    Shortened call to get_pao_gradient
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
!!    Included in pao_minimisation
!!   08:29, 2003/04/03 dave
!!    Changed to use pao_gradient for get_pao_gradient and get_electron_gradient
!!   09:14, 2003/04/10 dave
!!    Completely rewrote in a more transparent way (closely based on safemin in move_atoms.module)
!!   2008/05/25 ast
!!    Added timers
!!  SOURCE
!!
  subroutine line_minimise_pao( search_direction, fixed_potential, vary_mu, n_cg_L_iterations, &
       number_of_bands, tolerance, con_tolerance, mu, total_energy_0, expected_reduction, last_step, g_dot_sd )

    use datatypes
    use numbers
    use logicals
    use mult_module, ONLY: LNV_matrix_multiply, matM12, matM4
    use GenBlas, ONLY: copy, axpy, dot
    use SelfCon, ONLY: new_SC_potl
    use S_matrix_module, ONLY: get_S_matrix

    use GenComms, ONLY: gsum, my_barrier, cq_abort, inode, ionode
    ! Check on whether we need K found from L or whether we have it exactly
    use DiagModule, ONLY: diagon
    use global_module, ONLY: flag_vary_basis, ni_in_cell
    use primary_module, ONLY: bundle
    use support_spec_format, ONLY: supports_on_atom, coeff_array_size, coefficient_array, flag_paos_atoms_in_cell

    implicit none

    ! Passed variables

    logical :: vary_mu, fixed_potential, reset_L

    integer :: n_cg_L_iterations

    real(double) :: number_of_bands, tolerance, con_tolerance, mu
    real(double) :: total_energy_0
    real(double) :: expected_reduction, last_step, g_dot_sd
    real(double), dimension(:) :: search_direction

    ! Local variables
    integer :: lengthBlip, n_atoms

    real(double), save :: dE = 0.0_double ! Use this to guess initial step ?
    real(double), dimension(:), allocatable :: data_PAO0
    real(double), dimension(:), allocatable :: data_PAO
    real(double), dimension(:), allocatable :: data_full
    real(double) :: k0, k1, k2, k3, kmin, lambda
    real(double) :: e0, e1, e2, e3, electrons, tmp, energy_out, sum
    integer :: i,j, iter, part, memb, nsf1, npao1, iprim,l1,acz,m1

    logical :: done = .false. ! flag of line minimisation
    real(double), save :: kmin_last = 0.0_double

    write(io_lun,*) 'On entry to pao line_min, dE is ',dE, total_energy_0
    if(flag_paos_atoms_in_cell) then
       n_atoms = ni_in_cell
    else
       n_atoms = bundle%n_prim
    end if
    lengthBlip = coeff_array_size
    tmp = dot(lengthBlip,search_direction,1,search_direction,1)
    if(.NOT.flag_paos_atoms_in_cell) call gsum(tmp)
    !search_direction = search_direction/sqrt(tmp)
    write(io_lun,*) 'Searchdir: ',tmp
    !do i = 1,n_atoms
    !   do nsf1=1,nsf
    !      write(io_lun,*) 'Atom ',i,' supp ',nsf1,' coeffs ',supports_on_atom(i)%supp_func(nsf1)%coefficients
    !   end do
    !end do
    ! First, make a copy of the coefficients FOR THIS PRIMARY SET
    call start_timer(tmr_std_allocation)
    allocate(data_PAO0(lengthBlip))
    call stop_timer(tmr_std_allocation)
    data_PAO0 = coefficient_array
    ! We're assuming that we've ALREADY gone to a self-consistent ground state before arriving here
    iter = 1
    k1 = 0.0_double
    e1 = total_energy_0
    k2 = 0.0_double
    e2 = total_energy_0
    e3 = e2
    lambda = 2.0_double

    ! Loop to find a bracketing triplet
    if(dE== 0.0_double.or.kmin_last==0.0_double) then
      k3=8.0_double
    else
      !k3=0.5_double*dE/g_dot_sd
      k3=kmin_last
    endif

    done = .false.
    do while(.NOT.done)
       call copy( lengthBlip, data_PAO0, 1, coefficient_array, 1)
       call axpy( lengthBlip, k3, search_direction, 1, coefficient_array, 1 )
       ! Normalise
       do i = 1,n_atoms
          do nsf1=1,supports_on_atom(i)%nsuppfuncs
             sum = zero
             do npao1=1,supports_on_atom(i)%supp_func(nsf1)%ncoeffs
                sum = sum + supports_on_atom(i)%supp_func(nsf1)%coefficients(npao1)* &
                     supports_on_atom(i)%supp_func(nsf1)%coefficients(npao1)
             end do
             supports_on_atom(i)%supp_func(nsf1)%coefficients = supports_on_atom(i)%supp_func(nsf1)%coefficients/sqrt(sum)
          end do
       end do
       ! Find new self-consistent energy 
       ! 1. Get new S matrix (includes pao-to-grid transform)
       call get_S_matrix(inode, ionode)
       ! 2. If we're building K as 3LSL-2LSLSL, we need to make K now
       if(.NOT.diagon) call LNV_matrix_multiply(electrons, tmp, &
            doK, dontM1, dontM2, dontM3, dontM4, dontphi, dontE,0,0,0,0)
       reset_L = .true.
       ! 3. Get a new self-consistent potential and Hamiltonian
       ! I've not put a call to get_H_matrix here because it's currently in new_SC_potl
       call new_SC_potl( .false., con_tolerance, reset_L, fixed_potential, vary_mu, n_cg_L_iterations, &
            number_of_bands, tolerance, mu, e3)
       if(inode==ionode) write(io_lun,fmt='(2x,"In pao_min, iter ",i3," step and energy are ",2f15.10)') iter,k3,e3
       if(inode==ionode) write(io_lun,fmt='(" iter=", i3," k0,k1,k2,k3,kmin = ",5f15.8)') iter,k0,k1,k2,k3,kmin
       if (e3<e2) then ! We're still going down hill
          k1 = k2
          e1 = e2
          k2 = k3
          e2 = e3
          k3 = k3 * lambda  
          iter=iter+1
       else if(k2==0.0_double) then
          k3 = k3 / lambda
       else
          done = .true.
       endif
       if(k3 < very_small) call cq_abort('Step too small: line_minimise_pao failed!')
    end do
    ! Turn  basis variation back on
    ! Interpolate to find minimum.
    kmin = 0.5_double*(((k1*k1-k3*k3)*(e1-e2)-(k1*k1-k2*k2)*(e1-e3))/((k1-k3)*(e1-e2)-(k1-k2)*(e1-e3)))
    kmin_last = kmin
    if(inode==ionode) write(io_lun,*) 'In pao_min, bracketed - min from extrap: ',k1,k2,k3,kmin
    if(inode==ionode) write(io_lun,*) 'In pao_min, bracketed - energies: ',e1,e2,e3
    ! Change blips: start from blip0
    call copy( lengthBlip, data_PAO0, 1, coefficient_array, 1)
    call axpy( lengthBlip, kmin, search_direction, 1, coefficient_array, 1 )
    do i = 1,n_atoms
       do nsf1=1,supports_on_atom(i)%nsuppfuncs
          sum = zero
          do npao1=1,supports_on_atom(i)%supp_func(nsf1)%ncoeffs
             sum = sum + supports_on_atom(i)%supp_func(nsf1)%coefficients(npao1)* &
                  supports_on_atom(i)%supp_func(nsf1)%coefficients(npao1)
          end do
          supports_on_atom(i)%supp_func(nsf1)%coefficients = supports_on_atom(i)%supp_func(nsf1)%coefficients/sqrt(sum)
       end do
    end do
    ! Find new self-consistent energy 
    ! 1. Get new S matrix (includes pao-to-grid transform)
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
       call copy( lengthBlip, data_PAO0, 1, coefficient_array, 1)
       call axpy( lengthBlip, kmin, search_direction, 1, coefficient_array, 1 )
       do i = 1,n_atoms
          do nsf1=1,supports_on_atom(i)%nsuppfuncs
             sum = zero
             do npao1=1,supports_on_atom(i)%supp_func(nsf1)%ncoeffs
                sum = sum + supports_on_atom(i)%supp_func(nsf1)%coefficients(npao1)* &
                     supports_on_atom(i)%supp_func(nsf1)%coefficients(npao1)
             end do
             supports_on_atom(i)%supp_func(nsf1)%coefficients = supports_on_atom(i)%supp_func(nsf1)%coefficients/sqrt(sum)
          end do
       end do
       ! Find new self-consistent energy 
       ! 1. Get new S matrix (includes pao-to-grid transform)
       call get_S_matrix(inode, ionode)
       ! 2. If we're building K as 3LSL-2LSLSL, we need to make K now
       if(.NOT.diagon) call LNV_matrix_multiply(electrons, tmp, &
            doK, dontM1, dontM2, dontM3, dontM4, dontphi, dontE,0,0,0,0)
       reset_L = .true.
       ! 3. Get a new self-consistent potential and Hamiltonian
       call new_SC_potl( .false., con_tolerance, reset_L, fixed_potential, vary_mu, n_cg_L_iterations, &
            number_of_bands, tolerance, mu, energy_out)
    end if
    if(inode==ionode) write(io_lun,fmt='(2x,"In pao_min, at exit energy is ",2f15.10)') energy_out
    dE = total_energy_0 - energy_out
    write(io_lun,*) 'On exit from pao line_min, dE is ',dE,total_energy_0,energy_out
    total_energy_0 = energy_out
    deallocate(data_PAO0)
    return
  end subroutine line_minimise_pao
!!***

!!****f* pao_minimisation/build_PAO_coeff_grad *
!!
!!  NAME 
!!   build_PAO_coeff_grad
!!  USAGE
!! 
!!  PURPOSE
!!   Builds the gradient of energy wrt PAO coefficients
!!  INPUTS
!!   flag: if 1, do G.dS, if 2 do K.dH, if 3 do both.
!! 
!!   The gradient is calculated for atoms in this processor's primary set.
!!   So if we're storing coefficients for ALL atoms locally (shown by the
!!   flag_paos_atoms_in_cell), we need to do a gsum at the end (having stored 
!!   the coefficients for each primary atom in the appropriate part of the array)
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   2004 sometime
!!  MODIFICATION HISTORY
!!   2006/06/21 08:18 dave
!!    Changed for variable NSF and new basis storage scheme
!!   2008/05/25 ast
!!    Added timers
!!  SOURCE
!!

  subroutine build_PAO_coeff_grad(flag)

    use datatypes
    use logicals
    use numbers
    use DiagModule, ONLY: diagon
    use GenBlas, only: dot
    use primary_module, ONLY: bundle
    use cover_module, ONLY: BCS_parts
    use matrix_module, ONLY: matrix, matrix_halo
    use matrix_data, ONLY: mat,Srange,Hrange,halo,dSrange,dHrange
    use mult_module, ONLY: LNV_matrix_multiply, matM12, matM4, matS, matK, matdH, matdS, &
         return_matrix_value_pos, matrix_pos
    use GenComms, ONLY: myid, my_barrier, cq_abort,gsum
    use support_spec_format, ONLY: support_gradient, support_elec_gradient, flag_paos_atoms_in_cell,&
         grad_coeff_array, elec_grad_coeff_array, coeff_array_size
    use global_module, ONLY: iprint_minE

    ! Passed variables
    integer :: flag

    ! Local variables
    integer :: iprim, part, memb, nsf1, npao1, point, npao2, stat, i1, i2, &
         nab,gcspart, gcspartH,ist, nsfi, npaoi, S_nd_nabs,H_nd_nabs, count, atom_no, tmpc
    real(double) :: sum, e1, e2
    real(double), dimension(:),allocatable :: tmpS, tmpM12,tmpM4,tmpK,tmpdS,tmpdH

    if(.NOT.diagon) call LNV_matrix_multiply(e1, e2, &
         doK, doM1, doM2, dontM3, doM4, dontphi, dontE,matM12,0,matM4,0)
    ! We should have the elements built by H_matrix_module and S_matrix_module
    ! Now we take the sum over j\beta (nsf2 = \beta; neigh = j)
    iprim = 0
    call start_timer(tmr_std_matrices)
    do part = 1, bundle%groups_on_node
       if(bundle%nm_nodgroup(part) > 0) then
          do memb=1,bundle%nm_nodgroup(part) ! Select i
             nsfi = mat(part,Srange)%ndimi(memb)
             npaoi = mat(part,dSrange)%ndimi(memb)
             if(iprint_minE>2) write(io_lun,*) myid,' atom, nsf, npao: ',memb,nsfi,npaoi
             if(memb<bundle%nm_nodgroup(part)) then
                S_nd_nabs = mat(part,Srange)%i_nd_acc(memb+1)-mat(part,Srange)%i_nd_acc(memb)
                if(iprint_minE>2) write(io_lun,*) myid,' 1 nd: ',&
                     mat(part,Srange)%i_nd_acc(memb+1),mat(part,Srange)%i_nd_acc(memb)
             else 
                S_nd_nabs = mat(part,Srange)%part_nd_nabs-mat(part,Srange)%i_nd_acc(memb)+1
                if(iprint_minE>2) write(io_lun,*) myid,' 2 nd: ',&
                     mat(part,Srange)%part_nd_nabs,mat(part,Srange)%i_nd_acc(memb)
             end if
             if(iprint_minE>2) write(io_lun,*) 'S_nd_nabs: ', S_nd_nabs, mat(part,Srange)%n_nab(memb)
             S_nd_nabs = S_nd_nabs/nsfi
             if(iprint_minE>2) write(io_lun,*) 'S_nd_nabs: ', S_nd_nabs, mat(part,Srange)%n_nab(memb)
             call start_timer(tmr_std_allocation)
             allocate(tmpS(nsfi*S_nd_nabs))
             allocate(tmpM12(nsfi*S_nd_nabs))
             allocate(tmpM4(nsfi*S_nd_nabs))
             allocate(tmpdS(npaoi*S_nd_nabs))
             call stop_timer(tmr_std_allocation)
             if(memb<bundle%nm_nodgroup(part)) then
                H_nd_nabs = mat(part,Hrange)%i_nd_acc(memb+1)-mat(part,Hrange)%i_nd_acc(memb)
             else if(bundle%nm_nodgroup(part)>1) then
                H_nd_nabs = mat(part,Hrange)%part_nd_nabs-mat(part,Hrange)%i_nd_acc(memb)+1
             else
                H_nd_nabs = mat(part,Hrange)%part_nd_nabs
             end if
             if(iprint_minE>2) write(io_lun,*) myid,' H_nd_nabs: ', H_nd_nabs, mat(part,Hrange)%n_nab(memb)
             H_nd_nabs = H_nd_nabs/nsfi
             if(iprint_minE>2) write(io_lun,*) myid,' H_nd_nabs: ', H_nd_nabs, mat(part,Hrange)%n_nab(memb)
             call start_timer(tmr_std_allocation)
             allocate(tmpK(nsfi*H_nd_nabs))
             allocate(tmpdH(npaoi*H_nd_nabs))
             call stop_timer(tmr_std_allocation)
             tmpS = zero
             tmpM12 = zero
             tmpM4 = zero
             tmpdS = zero
             tmpK = zero
             tmpdH = zero
             iprim = iprim + 1
             if(flag_paos_atoms_in_cell) then
                atom_no = bundle%ig_prim(iprim)
             else
                atom_no = iprim
             end if
             count = 0
             tmpc = 0
             do nab = 1,mat(part,Srange)%n_nab(memb)
                ist = mat(part,Srange)%i_acc(memb)+nab-1
                gcspart = BCS_parts%icover_ibeg(mat(part,Srange)%i_part(ist))+mat(part,Srange)%i_seq(ist)-1
                do i2 = 1,mat(part,Srange)%ndimj(ist)
                   count = count + 1
                   do i1=1,nsfi
                      point = matrix_pos(matM12,iprim,halo(Srange)%i_halo(gcspart),i1,i2)
                      if((i1-1)*S_nd_nabs + count>nsfi*S_nd_nabs) call cq_abort("Overflow1 in buildPAOGrad ", &
                           (i1-1)*S_nd_nabs + count,nsfi*S_nd_nabs)
                      tmpM12((i1-1)*S_nd_nabs + count) = return_matrix_value_pos(matM12,point)
                      tmpM4((i1-1)*S_nd_nabs + count) = return_matrix_value_pos(matM4,point)
                      tmpS((i1-1)*S_nd_nabs + count) = return_matrix_value_pos(matS,point)
                      tmpc = tmpc+1
                   end do
                   do i1=1,npaoi
                      point = matrix_pos(matdS,iprim,halo(dSrange)%i_halo(gcspart),i1,i2)
                      if((i1-1)*S_nd_nabs + count>npaoi*S_nd_nabs) call cq_abort("Overflow2 in buildPAOGrad ", &
                           (i1-1)*S_nd_nabs + count,npaoi*S_nd_nabs)
                      tmpdS((i1-1)*S_nd_nabs + count) = return_matrix_value_pos(matdS,point)
                   end do
                end do
             end do
             if(iprint_minE>2) write(io_lun,*) 'tmpc is ',tmpc
             count = 0
             do nab = 1,mat(part,Hrange)%n_nab(memb)
                ist = mat(part,Hrange)%i_acc(memb)+nab-1
                gcspart = BCS_parts%icover_ibeg(mat(part,Hrange)%i_part(ist))+mat(part,Hrange)%i_seq(ist)-1
                do i2 = 1,mat(part,Hrange)%ndimj(ist)
                   count = count + 1
                   do i1=1,nsfi
                      point = matrix_pos(matK,iprim,halo(Hrange)%i_halo(gcspart),i1,i2)
                      tmpK((i1-1)*H_nd_nabs + count) = return_matrix_value_pos(matK,point)
                   end do
                   do i1=1,npaoi
                      point = matrix_pos(matdH,iprim,halo(dHrange)%i_halo(gcspart),i1,i2)
                      tmpdH((i1-1)*H_nd_nabs + count) = return_matrix_value_pos(matdH,point)
                   end do
                end do
             end do
             !write(myid+16,*) 'Atom: ',iprim,npaoi,nsfi
             do i1 = 1,nsfi
                do i2 = 1,npaoi
                   ! First do G.dS
                   if(iprint_minE>2) write(io_lun,*) 'atom, nsf, npao, grad: ',atom_no,i1, &
                        i2, support_gradient(atom_no)%supp_func(i1)%coefficients(i2)
                   if(flag==GdS.OR.flag==full) then
                      support_gradient(atom_no)%supp_func(i1)%coefficients(i2) =  &
                           support_gradient(atom_no)%supp_func(i1)%coefficients(i2)  - &
                           four*dot(S_nd_nabs,tmpM12((i1-1)*S_nd_nabs+1:),1,tmpdS((i2-1)*S_nd_nabs+1:),1)
                      !write(myid+16,*) 'dS: ',tmpdS((i2-1)*S_nd_nabs+1:i2*S_nd_nabs)
                      !write(myid+16,*) 'M12: ',tmpM12((i1-1)*S_nd_nabs+1:i1*S_nd_nabs)
                      !write(myid+16,*) 'grad: ',support_gradient(atom_no)%supp_func(i1)%coefficients(i2)
                      !write(myid+16,*) 'grad: ',four*dot(S_nd_nabs,tmpM12((i1-1)*S_nd_nabs+1:),1,tmpdS((i2-1)*S_nd_nabs+1:),1)
                   end if
                   if(iprint_minE>2) write(io_lun,*) 'atom, nsf, npao, grad: ',atom_no,i1, &
                        i2, support_gradient(atom_no)%supp_func(i1)%coefficients(i2)
                   ! Now do K.dH
                   if(flag==KdH.OR.flag==full) then
                      support_gradient(atom_no)%supp_func(i1)%coefficients(i2) =  &
                           support_gradient(atom_no)%supp_func(i1)%coefficients(i2)  - &
                           four*dot(H_nd_nabs,tmpK((i1-1)*H_nd_nabs+1:),1,tmpdH((i2-1)*H_nd_nabs+1:),1)
                      !write(myid+16,*) 'dH: ',tmpdH((i2-1)*H_nd_nabs+1:i2*H_nd_nabs)
                      !write(myid+16,*) 'K: ',tmpK((i1-1)*H_nd_nabs+1:i1*H_nd_nabs)
                      !write(myid+16,*) 'grad: ',four*dot(H_nd_nabs,tmpK((i1-1)*H_nd_nabs+1:),1,tmpdH((i2-1)*H_nd_nabs+1:),1)
                      !write(myid+16,*) 'grad: ',support_gradient(atom_no)%supp_func(i1)%coefficients(i2)
                   end if
                   if(iprint_minE>2) write(io_lun,*) 'atom, nsf, npao, grad: ',atom_no,i1, &
                        i2, support_gradient(atom_no)%supp_func(i1)%coefficients(i2)
                   ! Electron gradient
                   if(.NOT.diagon) then ! No problems with electron number when diagonalising
                      support_elec_gradient(atom_no)%supp_func(i1)%coefficients(i2) =  &
                           support_elec_gradient(atom_no)%supp_func(i1)%coefficients(i2)  - &
                           dot(S_nd_nabs,tmpM4((i1-1)*S_nd_nabs+1:),1,tmpdS((i2-1)*S_nd_nabs+1:),1)
                   end if
                end do ! i2=npao
             end do ! i1 = nsf
             call start_timer(tmr_std_allocation)
             deallocate(tmpS, tmpM12,tmpM4,tmpK,tmpdS,tmpdH)
             call stop_timer(tmr_std_allocation)
             !call cq_abort("Stopping now")
          end do ! memb=bundle%nm_nodgroup(part)
       end if
    end do ! part=bundle%groups_on_node    
    call stop_timer(tmr_std_matrices)
    if(flag_paos_atoms_in_cell) then
       call my_barrier
       call gsum(grad_coeff_array,coeff_array_size)
       call gsum(elec_grad_coeff_array,coeff_array_size)
    end if
    !write(15,*) grad_coeff_array
    return
  end subroutine build_PAO_coeff_grad

                   !*** Test gradients ***!
                   ! Shift coefficient
                   !if(inode==ionode) tmp = 0.001_double*supports_on_atom(bundle%ig_prim(iprim))%supp_func(nsf1)%coefficients(npao1)
                   ! supports_on_atom(bundle%ig_prim(iprim))%supp_func(nsf1)%coefficients(npao1) = &
                   !      supports_on_atom(bundle%ig_prim(iprim))%supp_func(nsf1)%coefficients(npao1) + tmp
                   !write(io_lun,*) 'On this proc, global(iprim) is ',iprim,bundle%ig_prim(iprim)
                   !%%!  tmp = 0.001_double*supports_on_atom(1)%supp_func(nsf1)%coefficients(npao1)
                   !%%!  supports_on_atom(1)%supp_func(nsf1)%coefficients(npao1) = &
                   !%%!       supports_on_atom(1)%supp_func(nsf1)%coefficients(npao1) + tmp
                   !%%!  ! Get gradient
                   !%%!  ! 1. Generate data_dS
                   !%%!  call get_S_matrix(support, inode, ionode, ntwof, SUPPORT_SIZE)
                   !%%!  ! 3. Generate data_dH
                   !%%!  call PAO_to_grid(inode-1,support,ntwof,SUPPORT_SIZE)
                   !%%!  call get_H_matrix( iprint_minE, .true., fixed_potential, &
                   !%%!       total_energy_test, electrons, ewald_energy, core_correction, &
                   !%%!       potential, density, pseudopotential, &
                   !%%!       support, workspace_support, inode, ionode, &
                   !%%!       N_GRID_MAX, ntwof, SUPPORT_SIZE)
                   !%%!  call FindMinDM(n_cg_L_iterations, number_of_bands, vary_mu, &
                   !%%!       L_tolerance, mu, inode, ionode, .false., .false.)
                   !%%!  call get_energy(E2, core_correction, ewald_energy)
                   !%%!  !call new_SC_potl( .false., con_tolerance, &
                   !%%!  !     reset_L, fixed_potential, vary_mu, n_cg_L_iterations, &
                   !%%!  !     number_of_bands, tolerance, mu,&
                   !%%!  !     E2, ewald_energy, core_correction,&
                   !%%!  !     potential, pseudopotential, density,&
                   !%%!  !     support, workspace_support, workspace2_support,&
                   !%%!  !     inode, ionode,&
                   !%%!  !     N_GRID_MAX, ntwof, SUPPORT_SIZE)
                   !%%!  gradient(npao1,nsf1,iprim) = zero
                   !%%!  point = mat(part,Srange)%offset+mat(part,Srange)%i_acc(memb)
                   !%%!  sum = dot(nsf*mat(part,Srange)%n_nab(memb),data_M12(nsf1,1:,point:),1,data_dS(npao1,1:,point:),1)
                   !%%!  gradient(npao1,nsf1,iprim) = gradient(npao1,nsf1,iprim) + four*sum
                   !%%!  !tmpgrad = four*sum
                   !%%!  point = mat(part,Hrange)%offset+mat(part,Hrange)%i_acc(memb)
                   !%%!  sum = dot(nsf*mat(part,Hrange)%n_nab(memb),data_K(nsf1,1:,point:),1,data_dH(npao1,1:,point:),1)
                   !%%!  gradient(npao1,nsf1,iprim) = gradient(npao1,nsf1,iprim) + four*sum
                   !%%!  !H2 = data_H(nsf1,npao,point)
                   !%%!  !H2a = data_H(nsf1,3,point)
                   !%%!  !E2 = nl_energy
                   !%%!  !BE2 = band_energy
                   !%%!  tmpgrad = gradient(npao1,nsf1,iprim)
                   !%%!  ! Shift coefficient
                   !%%!  !supports_on_atom(bundle%ig_prim(iprim))%supp_func(nsf1)%coefficients(npao1) = &
                   !%%!  !     supports_on_atom(bundle%ig_prim(iprim))%supp_func(nsf1)%coefficients(npao1) - tmp
                   !%%!  supports_on_atom(1)%supp_func(nsf1)%coefficients(npao1) = &
                   !%%!       supports_on_atom(1)%supp_func(nsf1)%coefficients(npao1) - tmp
                   !%%!  ! Get gradient
                   !%%!  ! 1. Generate data_dS
                   !%%!  call get_S_matrix(support, inode, ionode, ntwof, SUPPORT_SIZE)
                   !%%!  ! 3. Generate data_dH
                   !%%!  call PAO_to_grid(inode-1,support,ntwof,SUPPORT_SIZE)
                   !%%!  call get_H_matrix( iprint_minE, .true., fixed_potential, &
                   !%%!       total_energy_test, electrons, ewald_energy, core_correction, &
                   !%%!       potential, density, pseudopotential, &
                   !%%!       support, workspace_support, inode, ionode, &
                   !%%!       N_GRID_MAX, ntwof, SUPPORT_SIZE)
                   !%%!  call FindMinDM(n_cg_L_iterations, number_of_bands, vary_mu, &
                   !%%!       L_tolerance, mu, inode, ionode, .false., .false.)
                   !%%!  call get_energy(E1, core_correction, ewald_energy)
                   !%%!  !call new_SC_potl( .false., con_tolerance, &
                   !%%!  !     reset_L, fixed_potential, vary_mu, n_cg_L_iterations, &
                   !%%!  !     number_of_bands, tolerance, mu,&
                   !%%!  !     E1, ewald_energy, core_correction,&
                   !%%!  !     potential, pseudopotential, density,&
                   !%%!  !     support, workspace_support, workspace2_support,&
                   !%%!  !     inode, ionode,&
                   !%%!  !     N_GRID_MAX, ntwof, SUPPORT_SIZE)
                   !write(io_lun,*) 'Numerical, analytic grad: ',(E2-E1)/tmp, 0.5_double*(tmpgrad+gradient(npao1,nsf1,iprim))
                   ! BE1 = band_energy
                   !%%! H1 = data_H(nsf1,1,point)
                   !%%! H1a = data_H(nsf1,3,point)
                   !%%! E1 = nl_energy
                   !write(io_lun,*) 'Numerical, analytic grad: ',(BE2-BE1)/tmp, 0.5_double*(tmpgrad+four*sum)
                   !write(io_lun,*) 'Numerical, analytic grad: ',(E2-E1)/tmp, 0.5_double*(tmpgrad+gradient(npao1,nsf1,iprim))
                   !%%! write(io_lun,*) 'M Numerical, analytic: ',(H2-H1)/tmp,data_dH(npao1,1,point)
                   !%%! write(io_lun\,*) 'M Numerical, analytic: ',(H2a-H1a)/tmp,data_dH(npao1,3,point)
                   !sum = dot(nsf*mat(part,Srange)%n_nab(memb),data_K(nsf1,1:,point:),1,data_dC_NL(npao1,1,point),1)
                   !gradient(npao1,nsf1,iprim) = gradient(npao1,nsf1,iprim) + sum
                   !sum = dot(nsf*mat(part,Srange)%n_nab(memb),data_K(nsf1,1:,point:),1,data_dHloc(npao1,1,point),1)
                   !gradient(npao1,nsf1,iprim) = gradient(npao1,nsf1,iprim) + sum


end module pao_minimisation
