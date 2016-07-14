! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module blip_gradient
! ------------------------------------------------------------------------------
! Code area 6: Minimising energy
! ------------------------------------------------------------------------------

!!****h* Conquest/blip_gradient *
!!  NAME
!!   blip_gradient
!!  PURPOSE
!!   Contains all subroutines relating to finding the gradient of
!!   energy etc with respect to blips
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   08:25, 2003/04/03 dave
!!  MODIFICATION HISTORY
!!   10:09, 13/02/2006 drb 
!!    Removed all explicit references to data_ variables and rewrote
!!    in terms of new
!!    matrix routines
!!   2008/05/25 ast
!!    Added timers
!!   26/10/09 drb
!!    Various changes to allow selection of on-site analytics or not
!!   2012/04/02 17:15 dave
!!    Changes for analytic S and KE matrices
!!   2014/09/15 18:30 lat
!!    fixed call start/stop_timer to timer_module (not timer_stdlocks_module !)
!!  SOURCE
!!
module blip_gradient

  use timer_stdclocks_module, only: tmr_std_eminimisation,   &
                                    tmr_std_allocation,      &
                                    tmr_std_matrices
  use timer_module,           only: start_timer, stop_timer
                                    
  implicit none

  ! RCS tag for object file identification
  character(len=80), save, private :: &
       RCSid = "$Id$"
!!***

contains

! -----------------------------------------------------------
! Subroutine get_blip_gradient
! -----------------------------------------------------------

  !!****f* blip_gradient/get_blip_gradient *
  !!
  !!  NAME 
  !!   get_blip_gradient
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Evaluates the gradient of energy wrt blip coefficients
  !!  INPUTS
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   C.M.Goringe
  !!  CREATION DATE
  !!   19/10/95
  !!  MODIFICATION HISTORY
  !!   14/08/2000 TM
  !!    Added new act_on_vector
  !!   15/01/2001 TM
  !!    Added new blip_to_grad
  !!   15/05/2001 dave
  !!    Added ROBODoc header, indented, changed act_on_vector
  !!   16/05/2001 dave
  !!    Simplified get_onsite_KE call
  !!   17/05/2001 dave
  !!    Shortened blip_to_grad call and inverse calls
  !!   23/05/2001 dave
  !!    Shortened subroutine call
  !!   08:10, 31/07/2002 dave
  !!    Changed call to get_support_gradient not to pass data_M12 - it's
  !!    now used from matrix_data
  !!   08:26, 2003/04/03 dave
  !!    Included in blip_gradient
  !!   08:34, 2003/04/03 dave
  !!    Tidied use of GenBlas
  !!   2006/09/19 08:22 dave
  !!    Use temporary support fn on grid
  !!   2008/05/25 ast
  !!    Added timers
  !!   2009/10/27 07:12 dave
  !!    Added flag for on-site analytic elements
  !!   2011/11/29 L.Tong
  !!    Added spin polarisation
  !!   2012/03/22 L.Tong
  !!   - changed spin implementation
  !!   2012/04/26 16:15 dave
  !!    Changes for analytic S and KE matrices
  !!   2016/07/13 18:30 nakata
  !!    Renamed H_on_supportfns -> H_on_atomf
  !!  SOURCE
  !!
  subroutine get_blip_gradient(inode, ionode)

    use datatypes
    use numbers
    use dimens,                      only: grid_point_volume
    use primary_module,              only: bundle
    use matrix_data,                 only: mat, Hrange, Srange
    use mult_module,                 only: matK, matM12,                   &
                                           scale_matrix_value,             &
                                           return_matrix_block_pos
    use GenBlas,                     only: scal, copy
    use pseudopotential_data,        only: non_local
    use support_spec_format,         only: supports_on_atom,               &
                                           support_gradient, grad_coeff_array
    use species_module,              only: nsf_species
    use set_bucket_module,           only: rem_bucket
    use calc_matrix_elements_module, only: act_on_vectors_new
    use blip_grid_transform_module,  only: blip_to_grad_new,               &
                                           inverse_blip_transform_new,     &
                                           inverse_blip_to_grad_new
    use global_module,               only: WhichPulay, BothPulay, sf,      &
                                           flag_onsite_blip_ana,           &
                                           nspin, spin_factor,             &
                                           area_minE,                      &
                                           flag_analytic_blip_int, PhiPulay
    use functions_on_grid,           only: H_on_atomf, gridfunctions, &
                                           fn_on_grid,                     &
                                           allocate_temp_fn_on_grid,       &
                                           free_temp_fn_on_grid
    use memory_module,               only: reg_alloc_mem, type_dbl,        &
                                           reg_dealloc_mem
    use GenComms,                    only: cq_abort

    implicit none

    ! Passed variables
    integer :: inode, ionode

    ! local variables
    real(double) :: ke, ke2, tmp
    integer      :: i, direction, np, nn, n1, n2, this_nsf, tmp_fn, &
                    stat, spin
    real(double), dimension(:,:,:), allocatable :: this_data_K

    call start_timer(tmr_std_eminimisation)
    
    ! first, get the gradient of energy wrt the support functions.
    ! this does NOT include the change in energy due to the change in the
    ! KE matrix, but it DOES include the change in KE due to the change
    ! in S matrix.
    !
    ! a note on units... get_support_gradient and get_non_local_gradient
    ! both return values in units of 'per grid point'. The blip inverse
    ! transform assumes units of 'per au', so we rescale by multiplying
    ! by the grid point volume just before the transform. Similar arguments
    ! apply to the kinetic energy.
    !
    ! a note on signs... both the above mentioned routines give -Grad,
    ! ie the gradient downwards.
    !
    ! We have support and h_on_atomf in support and workspace_support,
    ! so use workspace2_support as working area...
    ! This returns the support_gradient in workspace_support
    WhichPulay = BothPulay
    ! Do on-site if we have this
    if ((.NOT.flag_analytic_blip_int).AND.flag_onsite_blip_ana) then
       i = 1
       this_nsf = nsf_species(bundle%species(i))

       allocate(this_data_K(this_nsf,this_nsf,nspin), STAT=stat)
       if (stat /= 0) &
            call cq_abort("get_blip_gradient: failed to allocate this_data_K", &
                          stat)
       call reg_alloc_mem(area_minE, this_nsf * this_nsf * nspin, type_dbl)
       this_data_K = zero

       call start_timer(tmr_std_matrices)
       
       do np = 1, bundle%groups_on_node
          do nn = 1, bundle%nm_nodgroup(np)
             if (nsf_species(bundle%species(i)) /= this_nsf) then
                deallocate(this_data_K)
                this_nsf = nsf_species(bundle%species(i))
                allocate(this_data_K(this_nsf,this_nsf,nspin))
             end if

             this_data_K = zero
             do spin = 1, nspin
                call return_matrix_block_pos(matM12(spin),              &
                                             mat(np,Srange)%onsite(nn), &
                                             this_data_K(:,:,spin),     &
                                             this_nsf * this_nsf)
                ! spin_factor correction for spin non-polarised
                ! calculation is done within get_onsite_S_gradient
                call get_onsite_S_gradient(supports_on_atom(i),   &
                                           this_data_K(:,:,spin), &
                                           support_gradient(i),   &
                                           this_nsf, bundle%species(i))
                do n1 = 1, this_nsf
                   do n2 = 1, this_nsf
                      call scale_matrix_value(matM12(spin), np, nn, i, 0, &
                                              n1, n2, zero, 1)
                   end do
                end do
             end do ! spin

             i = i + 1

          end do ! do nn = 1, bundle%nm_nodgroup(np)
       end do ! do np = 1, bundle%groups_on_node

       call stop_timer(tmr_std_matrices)

       deallocate(this_data_K)
       if (stat /= 0) &
            call cq_abort("get_blip_gradient: failed to deallocate this_data_K", &
                          stat)
       call reg_dealloc_mem(area_minE, this_nsf * this_nsf * nspin, type_dbl)
    end if
    if(flag_analytic_blip_int) WhichPulay = PhiPulay
    ! use H_on_atomf(1) as work storage for support gradient
    call get_support_gradient(H_on_atomf(1), inode, ionode)
    ! now we need to accumulate the 'type 1' gradient of the non-local
    ! energy (see notes, 3/1/97)
    if (non_local) call get_non_local_gradient(H_on_atomf(1), inode, ionode)

    ! now we need to transform this into the blip basis first, apply
    ! the scaling for grid size,
    call scal(gridfunctions(H_on_atomf(1))%size, grid_point_volume,&
              gridfunctions(H_on_atomf(1))%griddata, 1)

    ! and then 
    call inverse_blip_transform_new(inode-1, H_on_atomf(1), &
                                    support_gradient, bundle%n_prim)

    ! Kinetic Energy; first, the change in onsite T matrix elements,
    ! and as we do so, clear the diagonal blocks of data K
    if (flag_onsite_blip_ana.AND.(.NOT.flag_analytic_blip_int)) then
       i = 1
       this_nsf = nsf_species(bundle%species(i))

       call start_timer(tmr_std_allocation)
       allocate(this_data_K(this_nsf,this_nsf,nspin), STAT=stat)
       if (stat /= 0) &
            call cq_abort("get_blip_gradient: failed to allocate this_data_K",&
                          stat)
       call reg_alloc_mem(area_minE, this_nsf * this_nsf * nspin, type_dbl)
       call stop_timer (tmr_std_allocation)
       
       call start_timer (tmr_std_matrices)
       do np = 1, bundle%groups_on_node
          do nn = 1,bundle%nm_nodgroup(np)

             if (nsf_species(bundle%species(i)) /= this_nsf) then
                call start_timer(tmr_std_allocation)
                deallocate(this_data_K)
                call stop_timer(tmr_std_allocation)
                this_nsf = nsf_species(bundle%species(i))
                call start_timer(tmr_std_allocation)
                allocate(this_data_K(this_nsf,this_nsf,nspin))
                call stop_timer(tmr_std_allocation)
             end if
             
             this_data_K = zero
             do spin = 1, nspin
                call return_matrix_block_pos(matK(spin),                &
                                             mat(np,Hrange)%onsite(nn), &
                                             this_data_K(:,:,spin),      &
                                             this_nsf * this_nsf)                
                ! spin_factor correction for spin non-polarsed
                ! calculation is done inside get_onsite_KE_gradient
                call get_onsite_KE_gradient(supports_on_atom(i), &
                                            this_data_K(:,:,spin),   &
                                            support_gradient(i), &
                                            this_nsf, bundle%species(i))
                do n1 = 1, this_nsf
                   do n2 = 1, this_nsf
                      call scale_matrix_value(matK(spin), np, nn, i, &
                                              0, n1, n2, zero, 1)
                   end do ! n2
                end do ! n1
             end do ! spin

             i = i + 1

          end do ! nn
       end do ! np

       call stop_timer(tmr_std_matrices)

       call start_timer(tmr_std_allocation)
       deallocate(this_data_K, STAT=stat)
       if (stat /= 0) &
            call cq_abort("get_blip_gradient: failed to deallocate this_data_K", &
                          stat)
       call reg_dealloc_mem(area_minE, this_nsf * this_nsf * nspin, type_dbl)
       call stop_timer(tmr_std_allocation)

    end if ! (flag_onsite_blip_ana)

    ! Now, for the offsite part, done one the integration grid.
    ! to do the KE bit we need to blip_to_grad transform, act with K,
    ! and inverse_blip_grad_transform back (this routine ACCUMULATES onto
    ! data_dblip)
    !
    ! because we now do the onsite KE analytically, we need to remove
    ! the onsite parts of this term. This can be done by removing the
    ! onsite parts of the K matrix (because KE = Tr[KT], this is
    ! equivalent to removing the onsite parts of T, and somewhat easier!
    !
    ! first, for x: do the blip_grad transform into workspace_support 
    if(.NOT.flag_analytic_blip_int) then
       tmp_fn = allocate_temp_fn_on_grid(sf)
       do direction = 1, 3 
          call blip_to_grad_new(inode-1, direction, H_on_atomf(1))
          ! now act with K - result into workspace2_support
          gridfunctions(tmp_fn)%griddata = zero
          do spin = 1, nspin
             ! accumulate the contribution from matK(spin) to tmp_fn
             call act_on_vectors_new(inode-1, rem_bucket(3), matK(spin), &
                  tmp_fn, H_on_atomf(1))
          end do
          ! apply the appropriate scaling factor
          call scal(gridfunctions(tmp_fn)%size,        &
               - spin_factor * grid_point_volume, &
               gridfunctions(tmp_fn)%griddata, 1)

          call inverse_blip_to_grad_new(inode-1, direction, tmp_fn, &
               support_gradient, bundle%n_prim)
       end do
       call free_temp_fn_on_grid(tmp_fn)
    end if
    call stop_timer(tmr_std_eminimisation)

    return
  end subroutine get_blip_gradient
  !!***


  ! -----------------------------------------------------------
  ! Subroutine get_electron_gradient
  ! -----------------------------------------------------------
  
  !!****f* blip_gradient/get_electron_gradient *
  !!
  !!  NAME 
  !!   get_electron_gradient
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   This subroutine obtains the gradient of the electron
  !!   number with respect to the support functions on
  !!   the grid. The gradient with respect to
  !! |phi_i> at grid point n is given by
  !!   24 sum_j (LSLij  - LSLSLij) |phi_j(n)>
  !! 
  !!   M4 = 24 (LSL - LSLSL) 
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   E.H.Hernandez
  !!  CREATION DATE
  !!   01/03/95
  !!  MODIFICATION HISTORY
  !!   12/9/95 EH to conform to new data storage
  !!   8/5/96 CG this routine becomes trivial...
  !!   7/7/98 DRB to use act_on_vector
  !!   19/9/2000 TM
  !!    Added new act_on_vector calls
  !!   15/05/2001 dave
  !!    Converted to F90, added ROBODoc and updated act_on_vector call
  !!   23/05/2001 dave
  !!    Shortened subroutine call
  !!   24/05/2001 dave
  !!    Bug fix - added use maxima_module
  !!   08:27, 2003/04/03 dave
  !!    Included in blip_gradient
  !!   08:34, 2003/04/03 dave
  !!    Tidied use of GenBlas
  !!   2011/12/05 L.Tong
  !!    Added spin polarisation
  !!   2012/03/23 L.Tong
  !!   - Changed spin implementation
  !!  SOURCE
  !!
  subroutine get_electron_gradient(support, electron_gradient, inode, &
                                   ionode)

    use datatypes
    use numbers
    use GenBlas,                     only: scal
    use dimens,                      only: grid_point_volume
    use set_bucket_module,           only: rem_bucket
    use calc_matrix_elements_module, only: act_on_vectors_new
    use blip_grid_transform_module,  only: inverse_blip_transform_new
    use mult_module,                 only: matM4
    use functions_on_grid,           only: gridfunctions, fn_on_grid
    use primary_module,              only: bundle
    use support_spec_format,         only: support_elec_gradient
    use global_module,               only: nspin, spin_factor

    implicit none

    ! Passed variables
    integer :: inode, ionode, spin
    integer :: support, electron_gradient

    gridfunctions(electron_gradient)%griddata = zero
    do spin = 1, nspin
       ! accumulate the contribution from matM4(spin) to
       ! electron_gradient
       call act_on_vectors_new(inode-1, rem_bucket(1), matM4(spin), &
                               electron_gradient, support)
    end do
    ! scale correctly
    call scal(gridfunctions(electron_gradient)%size, &
              - spin_factor * grid_point_volume,     &
              gridfunctions(electron_gradient)%griddata, 1)
    
    call inverse_blip_transform_new(inode-1, electron_gradient, &
                                    support_elec_gradient,      &
                                    bundle%n_prim)
    
    return
  end subroutine get_electron_gradient
  !!***

  ! -----------------------------------------------------------
  ! Subroutine get_support_gradient
  ! -----------------------------------------------------------
  
  !!****f* blip_gradient/get_support_gradient *
  !!
  !!  NAME 
  !!   get_support_gradient
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   This subroutine obtains the gradient of the total
  !!   energy with respect to the support functions on
  !!   the grid. The gradient of E with respect to
  !!   |phi_i> at grid point n is given by
  !!
  !!    dE/d|phi_i(n)> = 4 sum_j 
  !!            (      Kij H +          <- Type I
  !!                   M1ij - M2ij    )   <- Type II
  !!   |phi_j(n)>
  !! 
  !!    where H is the Hamiltonian operator (not its 
  !!    matrix representation), which must be made to 
  !!    act on |phi_j>, and the matrices M1 and M2 are
  !!    defined as
  !!    M1 = 3 LHL 
  !!    M2 = 2 LSLHL + LHLSL
  !!    respectively.
  !!
  !!    Notes 4/6/94, section 11.2 (pp24-34)
  !!    The Hamiltonian matrix here is the complete matrix, but the
  !!    Hamiltonian operator is (H_local - T); the non-local type I
  !!    derivative and the kinetic type I derivative have to be done
  !!    seperately.
  !!    See Notes 24/2/96 pp9-10 and Conquest paper IV for KE discussion,
  !!        Notes 31/1/97 pp1-5 for Non-Local energy
  !! 
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   E.H.Hernandez
  !!  CREATION DATE
  !!   01/03/95
  !!  MODIFICATION HISTORY
  !!   12/09/95 EH to conform to new data storage
  !!   02/05/96 EH to speed up and take advantage of
  !!             block data distribution.
  !!   03/05/96 CG to simplify above changes
  !!   19/9/2000 TM
  !!    Added new act_on_vector calls
  !!   15/05/2001 dave
  !!    F90, ROBODoc and act_on_vector sorted
  !!   21/05/2001 dave
  !!    Shortened variables passed
  !!   31/07/2002 dave
  !!    Changed to use data_M12 from matrix_data instead of passing it
  !!   08:27, 2003/04/03 dave
  !!    Included in blip_gradient
  !!   08:34, 2003/04/03 dave
  !!    Tidied use of GenBlas
  !!   2011/09/19 L.Tong
  !!    - Added spin polarisation
  !!    - Relies on both H_on_supportfns and H_dn_on_supportfns containing
  !!      the correct values upon entry to the subroutine.
  !!    - Upon exit, H_on_supportfns stores the TOTAL gradient.
  !!                 H_dn_on_supportfns stores the gradient for spin down component
  !!   2012/03/22 L.Tong
  !!   - Major rewrite
  !!   - Added support_gradient, so now the total support gradient is
  !!     stored in there
  !!   - H_on_support(1) and H_on_support(nspin) are now NOT modified
  !!     and used for storage of support gradient
  !!   - tmp_fn is now an array of dimension nspin, and the contribution
  !!     of type I and II gradients are accumulated using the fact that
  !!     act_on_vectors_new is accumulative. This improves code
  !!     efficiency.
  !!   2016/07/13 18:30 nakata
  !!    Renamed H_on_supportfns -> H_on_atomf
  !!  SOURCE
  !!
  subroutine get_support_gradient(support_gradient, inode, ionode)

    use datatypes
    use numbers
    use mult_module,                 only: matM12, matK
    use GenBlas,                     only: axpy, copy, scal
    use set_bucket_module,           only: rem_bucket
    use calc_matrix_elements_module, only: act_on_vectors_new
    use global_module,               only: flag_basis_set, PAOs,      &
                                           blips, WhichPulay,         &
                                           PhiPulay, BothPulay,       &
                                           SPulay, sf, nspin,         &
                                           spin_factor
    use functions_on_grid,           only: gridfunctions, fn_on_grid, &
                                           supportfns,                &
                                           H_on_atomf,           &
                                           allocate_temp_fn_on_grid,  &
                                           free_temp_fn_on_grid

    implicit none

    ! Shared variables
    integer :: support_gradient, inode, ionode

    ! Local variables
    integer :: spin
    integer, dimension(nspin) :: tmp_fn

    do spin = 1, nspin
       tmp_fn(spin) = allocate_temp_fn_on_grid(sf)
       gridfunctions(tmp_fn(spin))%griddata = zero
    end do

    ! 2011/09/20 L.Tong
    ! This subroutine is called on the premise that both H_on_atomf
    ! and H_dn_on_supportfns (for spin polarised calculations) alreay contain
    ! the correct values upon entry.

    ! first, act on h_on_atomf (which holds H |phi> on entry)
    ! with the K matrix) to get type I variation
    
    ! L.Tong: act_on_vectors_new is accumulative, we use this to our advantage
    if (WhichPulay == PhiPulay .or. WhichPulay == BothPulay) then
       do spin = 1, nspin
          call act_on_vectors_new(inode-1, rem_bucket(3), matK(spin), &
                                  tmp_fn(spin), H_on_atomf(spin))
       end do
    end if
    ! now act with M12 on support to get type II variation
    ! (accumulate)
    if (flag_basis_set == blips .and. (WhichPulay == SPulay .or. &
         WhichPulay == BothPulay)) then
       do spin = 1, nspin
          call act_on_vectors_new(inode-1, rem_bucket(1), &
                                  matM12(spin), tmp_fn(spin), &
                                  supportfns)
       end do
    end if
    ! multiply by the correct factor for each spin channel
    do spin = 1, nspin
       call scal(gridfunctions(tmp_fn(spin))%size, minus_two, &
                 gridfunctions(tmp_fn(spin))%griddata, 1)
    end do
    
    ! get the total support gradient
    gridfunctions(support_gradient)%griddata = zero
    do spin = 1, nspin
       call axpy(gridfunctions(support_gradient)%size, spin_factor, &
                 gridfunctions(tmp_fn(spin))%griddata, 1,           &
                 gridfunctions(support_gradient)%griddata, 1)
    end do
    ! free tmp_fn
    do spin = nspin, 1, -1
       call free_temp_fn_on_grid(tmp_fn(spin))
    end do

    return
  end subroutine get_support_gradient
  !!***


  ! -----------------------------------------------------------
  ! Subroutine get_non_local_gradient
  ! -----------------------------------------------------------
  
  !!****f* blip_gradient/get_non_local_gradient *
  !!
  !!  NAME 
  !!   get_non_local_gradient
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   accumulates the gradient of energy due to type 1
  !!   variations of the non-local energy onto the grid
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   C.M.Goringe
  !!  CREATION DATE
  !!   24/01/97
  !!  MODIFICATION HISTORY
  !!   07/07/98 by Dave Bowler to use act_on_vector
  !!   15/05/2001 dave
  !!    Changed act_on_vector, added ROBODoc header, indented
  !!   23/05/2001 dave
  !!    Shortened subroutine call
  !!   15:42, 04/02/2003 drb 
  !!    Changed to use temp_TCS (problem with transposes if NSF/=NCF)
  !!   08:15, 2003/04/04 dave
  !!    Included into blip_gradient
  !!   2011/11/29 L.Tong
  !!    Added spin polarisation
  !!    - Upon exit, support_gradient stores the accumulative TOTAL gradient
  !!      contributions from BOTH spin components
  !!   2011/12/13 07:53 dave
  !!    Added analytic evaluation of blip integrals
  !!   2012/03/20 L.Tong
  !!   - Rewrite spin implementation
  !!   - for spin non-polarised calculation, matU(1) now only stores
  !!     the correct values for spin up component, so it has a factor
  !!     of two less than that of the original implementation. This
  !!     leads to the need for corrections involving spin_factor.
  !!   2012/03/22 L.Tong
  !!   - Improve efficiency by absorbing spin_factor into matU, so
  !!     that factor corrections for spin non-polarised calculations
  !!     are not required.
  !!  SOURCE
  !!
  subroutine get_non_local_gradient(support_grad_grid, inode, ionode)

    use datatypes
    use numbers
    use mult_module,                 only: H_SP_SP, SP_trans, mult, &
                                           matrix_product,          &
                                           matrix_scale,            &
                                           matrix_transpose, matSC, &
                                           matCS, matU, matK,       &
                                           return_matrix_block_pos, &
                                           matrix_pos
    use set_bucket_module,           only: rem_bucket
    use calc_matrix_elements_module, only: act_on_vectors_new
    use functions_on_grid,           only: pseudofns
    use nlpf2blip,                   only: get_blipP, nlpf_on_atom
    use primary_module ,             only: bundle
    use support_spec_format,         only: supports_on_atom,        &
                                           support_gradient
    use group_module,                only: parts
    use primary_module,              only: bundle
    use cover_module,                only: BCS_parts
    use matrix_data,                 only: mat, SPrange, halo
    use species_module,              only: nsf_species, nlpf_species
    use global_module,               only: nspin, spin_factor,      &
                                           id_glob, species_glob,   &
                                           flag_analytic_blip_int

    implicit none

    ! Passed variables
    integer :: inode, ionode
    integer :: support_grad_grid

    ! Local variables
    integer      :: iprim, np, nn, spec, nab, ist, gcspart,       &
                    neigh_global_part, neigh_global_num,          &
                    neigh_species, wheremat, this_nsf, this_nlpf, &
                    j_in_halo, spin
    real(double) :: dx, dy, dz
    real(double), dimension(:,:), allocatable :: dataU

    ! First of all, find U (=K.SC)
    call matrix_transpose(matSC,matCS)
    do spin = 1, nspin
       call matrix_product(matK(spin), matCS, matU(spin), mult(H_SP_SP))
       call matrix_scale(minus_two * spin_factor, matU(spin))
    end do

    ! Finally, act on the non-local projectors with the scaled U
    if (flag_analytic_blip_int) then
       iprim = 0
       do np = 1, bundle%groups_on_node
          do nn = 1, bundle%nm_nodgroup(np)
             iprim = iprim + 1
             spec = bundle%species(iprim)
             this_nsf = nsf_species(spec)
             do nab = 1, mat(np,SPrange)%n_nab(nn) ! Loop over neighbours of atom
                ist = mat(np,SPrange)%i_acc(nn) + nab - 1
                ! Build the distances between atoms - needed for phases 
                gcspart = &
                     BCS_parts%icover_ibeg(mat(np,SPrange)%i_part(ist)) + &
                     mat(np,SPrange)%i_seq(ist) - 1
                ! Displacement vector
                dx = BCS_parts%xcover(gcspart) - bundle%xprim(iprim)
                dy = BCS_parts%ycover(gcspart) - bundle%yprim(iprim)
                dz = BCS_parts%zcover(gcspart) - bundle%zprim(iprim)
                ! We need to know the species of neighbour
                neigh_global_part = &
                     BCS_parts%lab_cell(mat(np,SPrange)%i_part(ist)) 
                neigh_global_num  = &
                     id_glob(parts%icell_beg(neigh_global_part) + &
                             mat(np,SPrange)%i_seq(ist) - 1)
                neigh_species = species_glob(neigh_global_num)
                this_nlpf = nlpf_species(neigh_species)
                j_in_halo = halo(SPrange)%i_halo(gcspart)
                wheremat = matrix_pos(matSC, iprim, j_in_halo, 1, 1)
                allocate(dataU(this_nsf,this_nlpf))
                do spin = 1, nspin
                   dataU = zero
                   ! note matU(1) for spin non-polarised version
                   ! already has factor of two (spin_factor)
                   ! multiplied to it. So no corrections needed here
                   call return_matrix_block_pos(matU(spin), wheremat, dataU, &
                                                this_nsf * this_nlpf)
                   ! get_blipP accumulates
                   call get_blipP(support_gradient(iprim),            &
                                  nlpf_on_atom(neigh_species), dataU, &
                                  iprim, j_in_halo, dx, dy, dz, spec, &
                                  neigh_species)
                end do
                deallocate(dataU)
             end do ! nab
          end do ! nn
       end do ! np
    else
       do spin = 1, nspin
          ! accmulate matU(spin).pseudofns
          call act_on_vectors_new(inode-1, rem_bucket(2), matU(spin), &
                                  support_grad_grid, pseudofns)
          ! note matU(1) for spin non-polarised version
          ! already has factor of two (spin_factor)
          ! multiplied to it. So no corrections needed here
       end do
    end if
    return
  end subroutine get_non_local_gradient
  !!***
  

  !!****f* blip_gradient/get_onsite_KE_gradient *
  !!
  !!  NAME 
  !!   get_onsite_KE_gradient
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   This routine accumulates onto the blip gradient the derivative in
  !!   KE due to the change in T wrt the blip coefficients
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   C.M.Goringe
  !!  CREATION DATE
  !!   18/11/96 
  !!  MODIFICATION HISTORY
  !!   16/05/2001 dave
  !!    F90, ROBODoc and:
  !!    REMOVED DEPENDENCE ON NSF=4
  !!   16/05/2001 dave
  !!    Simplified call
  !!   08:17, 2003/04/04 dave
  !!    Included into blip_gradient
  !!   2008/05/25 ast
  !!    Added timers
  !!   2012/02/20 L.Tong
  !!    Added spin polarisation
  !!   2012/03/22 L.Tong
  !!   - Changed spin implementation
  !!   - Changed this_data_K to a 2D array.
  !!   2012/04/26 17:04 dave
  !!    Small bug fix in indexing of this_data_K
  !!  SOURCE
  !!
  subroutine get_onsite_KE_gradient(this_data_blip, this_data_K, &
                                    this_blip_grad, this_nsf, spec)

    use datatypes
    use numbers
    use GenBlas,             only: axpy, copy, scal
    use blip,                only: blip_info
    use support_spec_format, only: support_function
    use GenComms,            only: cq_abort
    use global_module,       only: area_minE, spin_factor
    use memory_module,       only: reg_alloc_mem, type_dbl, reg_dealloc_mem

    implicit none

    ! Passed Variables
    integer                :: this_nsf, spec
    type(support_function) :: this_data_blip, this_blip_grad
    ! N.B. Here NSF is for ONE site only so is OK
    real(double), dimension(this_nsf,this_nsf) :: this_data_K

    ! Local Variables
    integer, parameter :: MAX_D=3
    integer      :: dx, dy, dz, offset, l, at, nsf1, nsf2, stat
    real(double) ::  FAC(0:MAX_D), D2FAC(0:MAX_D)
    real(double), allocatable, dimension(:) :: work1, work2, work3, &
                                               work4, work5, work6

    call start_timer(tmr_std_allocation)
    allocate(work1(blip_info(spec)%FullArraySize*this_nsf), &
         work2(blip_info(spec)%FullArraySize*this_nsf), &
         work3(blip_info(spec)%FullArraySize*this_nsf), &
         work4(blip_info(spec)%FullArraySize*this_nsf), &
         work5(blip_info(spec)%FullArraySize*this_nsf), &
         work6(blip_info(spec)%FullArraySize*this_nsf), &
         STAT=stat)
    if (stat /= 0) call cq_abort("Error allocating arrays for onsite &
         &KE blip grad: ", blip_info(spec)%FullArraySize,this_nsf)
    call reg_alloc_mem (area_minE, 6 * blip_info(spec)%FullArraySize *&
         this_nsf, type_dbl)
    call stop_timer (tmr_std_allocation)
    ! first, we copy the blip functions for this atom onto a cubic grid;
    ! we make this grid 'too big' in order to have a fast routine below.
    !
    ! To get gradients, we start by acting with the K matrix upon the
    ! coefficients

    FAC(0) = 151.0_double/140.0_double
    FAC(1) = 1191.0_double/2240.0_double
    FAC(2) = 3.0_double/56.0_double
    FAC(3) = 1.0_double/2240.0_double
    D2FAC(0) = 3.0_double/2.0_double
    D2FAC(1) = -9.0_double/32.0_double
    D2FAC(2) = -18.0_double/40.0_double
    D2FAC(3) = -3.0_double/160.0_double

    work1 = zero
    offset = blip_info(spec)%BlipArraySize+1

    do dx = -blip_info(spec)%BlipArraySize, blip_info(spec)%BlipArraySize
       do dy = -blip_info(spec)%BlipArraySize, blip_info(spec)%BlipArraySize
          do dz = -blip_info(spec)%BlipArraySize, blip_info(spec)%BlipArraySize
             l = blip_info(spec)%blip_number(dx,dy,dz)
             if (l .ne. 0) then
                at = (((dz + offset) * blip_info(spec)%OneArraySize + &
                     (dy + offset)) * blip_info(spec)%OneArraySize + &
                     (dx + offset)) * this_nsf
                do nsf1 = 1,this_nsf
                   work1(nsf1 + at) = zero
                   do nsf2 = 1, this_nsf
                      work1(nsf1 + at) = work1(nsf1 + at) + &
                           this_data_K(nsf2,nsf1) *&
                           this_data_blip%supp_func(nsf2)%coefficients(l)
                   enddo
                enddo
             end if
          end do
       end do
    end do

    ! now, for each direction in turn, we need to apply a number of
    ! 'spreading operations'. Do z first... put blip(z) in 2, and
    ! del2blip(z) in 3

    call copy (blip_info(spec)%FullArraySize * this_nsf, work1, 1, &
         work2, 1)
    call scal (blip_info(spec)%FullArraySize * this_nsf, FAC(0), &
         work2, 1)
    do dz = 1, MAX_D
       offset = dz * blip_info(spec)%OneArraySize * &
            blip_info(spec)%OneArraySize * this_nsf
       call axpy ((blip_info(spec)%FullArraySize * this_nsf - offset),&
            FAC(dz), work1(1:), 1, work2(1 + offset:), 1 )
       call axpy ((blip_info(spec)%FullArraySize * this_nsf - offset),&
            FAC(dz), work1(1 + offset:), 1, work2(1:), 1 )
    end do

    call copy (blip_info(spec)%FullArraySize * this_nsf, work1, 1, &
         work3, 1)
    call scal (blip_info(spec)%FullArraySize * this_nsf, D2FAC(0), &
         work3, 1)
    do dz = 1, MAX_D
       offset = dz * blip_info(spec)%OneArraySize * &
            blip_info(spec)%OneArraySize * this_nsf
       call axpy ((blip_info(spec)%FullArraySize * this_nsf - offset),&
            D2FAC(dz), work1(1:), 1, work3(1 + offset:), 1)
       call axpy ((blip_info(spec)%FullArraySize * this_nsf - offset),&
            D2FAC(dz), work1(1 + offset:), 1, work3(1:), 1)
    end do

    ! now do y : put blip(y).blip(z) in 4,
    ! blip(y).del2blip(z) + del2blip(y).blip(z)  in 5

    call copy (blip_info(spec)%FullArraySize * this_nsf, work2, 1, &
         work4, 1)
    call scal (blip_info(spec)%FullArraySize * this_nsf, FAC(0), &
         work4, 1)
    do dy = 1, MAX_D
       offset = dy * blip_info(spec)%OneArraySize * this_nsf
       call axpy ((blip_info(spec)%FullArraySize * this_nsf - offset),&
            FAC(dy), work2(1:), 1, work4(1 + offset:), 1 )
       call axpy ((blip_info(spec)%FullArraySize * this_nsf - offset),&
            FAC(dy), work2(1 + offset:), 1, work4(1:), 1 )
    end do

    call copy (blip_info(spec)%FullArraySize * this_nsf, work2, 1, &
         work5, 1)
    call scal (blip_info(spec)%FullArraySize * this_nsf, D2FAC(0), &
         work5, 1)
    do dy = 1, MAX_D
       offset = dy * blip_info(spec)%OneArraySize * this_nsf
       call axpy ((blip_info(spec)%FullArraySize * this_nsf - offset),&
            D2FAC(dy), work2(1:), 1, work5(1 + offset:), 1 )
       call axpy ((blip_info(spec)%FullArraySize * this_nsf-offset), &
            D2FAC(dy), work2(1 + offset:), 1, work5(1:), 1 )
    end do

    call axpy (blip_info(spec)%FullArraySize * this_nsf, FAC(0), &
         work3, 1, work5, 1)
    do dy = 1, MAX_D
       offset = dy * blip_info(spec)%OneArraySize * this_nsf
       call axpy ((blip_info(spec)%FullArraySize * this_nsf - offset),&
            FAC(dy), work3(1:), 1, work5(1 + offset:), 1)
       call axpy ((blip_info(spec)%FullArraySize * this_nsf - offset),&
            FAC(dy), work3(1 + offset:), 1, work5(1:), 1)
    end do

    ! and x - put it all into 6

    call copy (blip_info(spec)%FullArraySize * this_nsf, work5, 1, &
         work6, 1)
    call scal (blip_info(spec)%FullArraySize * this_nsf, FAC(0), &
         work6, 1)
    do dx = 1, MAX_D
       offset = dx * this_nsf
       call axpy ((blip_info(spec)%FullArraySize * this_nsf - offset),&
            FAC(dx), work5(1:), 1, work6(1 + offset:), 1)
       call axpy ((blip_info(spec)%FullArraySize * this_nsf - offset),&
            FAC(dx), work5(1 + offset:), 1, work6(1:), 1)
    end do

    call axpy (blip_info(spec)%FullArraySize * this_nsf, D2FAC(0), &
         work4, 1, work6, 1)
    do dx = 1, MAX_D
       offset = dx * this_nsf
       call axpy ((blip_info(spec)%FullArraySize * this_nsf - offset),&
            D2FAC(dx), work4(1:), 1, work6(1 + offset:), 1)
       call axpy ((blip_info(spec)%FullArraySize * this_nsf - offset),&
            D2FAC(dx), work4(1 + offset:), 1, work6(1:), 1)
    end do

    ! now accumulate work6 onto the gradient
    call scal(blip_info(spec)%FullArraySize * this_nsf, spin_factor * &
              blip_info(spec)%SupportGridSpacing, work6, 1)

    offset = blip_info(spec)%BlipArraySize + 1
    do dx = - blip_info(spec)%BlipArraySize, blip_info(spec)%BlipArraySize
       do dy = - blip_info(spec)%BlipArraySize, blip_info(spec)%BlipArraySize
          do dz = - blip_info(spec)%BlipArraySize, blip_info(spec)%BlipArraySize
             l = blip_info(spec)%blip_number(dx, dy, dz)
             if (l .ne. 0) then
                at = (((dz + offset) * blip_info(spec)%OneArraySize + &
                     (dy + offset)) * blip_info(spec)%OneArraySize + &
                     (dx + offset)) * this_nsf
                do nsf1 = 1, this_nsf
                   this_blip_grad%supp_func(nsf1)%coefficients(l) = &
                        this_blip_grad%supp_func(nsf1)%coefficients(l) - &
                        work6(nsf1 + at)
                enddo
             end if
          end do
       end do
    end do

    call start_timer (tmr_std_allocation)
    deallocate (work1, work2, work3, work4, work5, work6, STAT=stat)
    if (stat /= 0) &
         call cq_abort ("Error deallocating arrays for onsite KE blip grad: ", &
         blip_info(spec)%FullArraySize,this_nsf)
    call reg_dealloc_mem (area_minE, 6 * blip_info(spec)%FullArraySize * &
         this_nsf, type_dbl)
    call stop_timer (tmr_std_allocation)
    return
  end subroutine get_onsite_KE_gradient
!!***


  !!****f* blip_gradient/get_onsite_S_gradient *
  !!
  !!  NAME 
  !!   get_onsite_S_gradient
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   This routine accumulates onto the blip gradient the derivative in
  !!   energy due to the change in onsite S wrt the blip coefficients
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   26/10/09
  !!  MODIFICATION HISTORY
  !!   2011/11/15 08:04 dave
  !!    Changes to blip data
  !!   2012/02/20 L.Tong
  !!    Added spin polarisation
  !!   2012/03/22 L.Tong
  !!   - Changed spin implementation
  !!   2012/04/26 17:04 dave
  !!    Small bug fix in indexing of this_data_M
  !!  SOURCE
  !!
  subroutine get_onsite_S_gradient(this_data_blip, this_data_M, &
                                   this_blip_grad, this_nsf, spec)

    use datatypes
    use numbers
    use GenBlas,             only: axpy, copy, scal
    use blip,                only: blip_info
    use support_spec_format, only: support_function
    use GenComms,            only: cq_abort
    use global_module,       only: area_minE, spin_factor
    use memory_module,       only: reg_alloc_mem, type_dbl, reg_dealloc_mem

    implicit none

    ! Passed Variables
    integer :: this_nsf, spec
    ! N.B. Here NSF is for ONE site only so is OK
    real(double) :: this_data_M(this_nsf,this_nsf)
    type(support_function) :: this_data_blip, this_blip_grad

    ! Local Variables
    integer, parameter :: MAX_D=3
    real(double) ::  FAC(0:MAX_D)
    real(double), allocatable, dimension(:) :: work1, work2, work4, work6
    integer :: dx, dy, dz, offset, l, at, nsf1, nsf2, stat
    real(double) :: sum

    call start_timer (tmr_std_allocation)
    allocate (work1(blip_info(spec)%FullArraySize*this_nsf),&
         work2(blip_info(spec)%FullArraySize*this_nsf),&
         work4(blip_info(spec)%FullArraySize*this_nsf),&
         work6(blip_info(spec)%FullArraySize*this_nsf), STAT=stat)
    work1 = zero
    work2 = zero
    work4 = zero
    work6 = zero
    if (stat /= 0) &
         call cq_abort ("Error allocating arrays for onsite KE blip grad: ", &
         blip_info(spec)%FullArraySize, this_nsf)
    call reg_alloc_mem (area_minE, 4 * blip_info(spec)%FullArraySize * &
         this_nsf, type_dbl)
    call stop_timer (tmr_std_allocation)

    ! first, we copy the blip functions for this atom onto a cubic grid;
    ! we make this grid 'too big' in order to have a fast routine below.
    !
    ! To get gradients, we start by acting with the K matrix upon the
    ! coefficients

    FAC(0) = 151.0_double / 140.0_double
    FAC(1) = 1191.0_double / 2240.0_double
    FAC(2) = 3.0_double / 56.0_double
    FAC(3) = 1.0_double / 2240.0_double

    work1 = zero
    offset = blip_info(spec)%BlipArraySize+1

    do dx = - blip_info(spec)%BlipArraySize, blip_info(spec)%BlipArraySize
       do dy = - blip_info(spec)%BlipArraySize, blip_info(spec)%BlipArraySize
          do dz = - blip_info(spec)%BlipArraySize, blip_info(spec)%BlipArraySize
             l = blip_info(spec)%blip_number(dx,dy,dz)
             if (l .ne. 0) then
                at = (((dz + offset) * blip_info(spec)%OneArraySize + &
                     (dy + offset)) * blip_info(spec)%OneArraySize + &
                     (dx + offset)) * this_nsf
                do nsf1 = 1, this_nsf
                   work1(nsf1 + at) = zero
                   do nsf2 = 1, this_nsf
                      work1(nsf1 + at) = work1(nsf1 + at) + &
                           this_data_M(nsf2,nsf1) * &
                           this_data_blip%supp_func(nsf2)%coefficients(l)
                   enddo
                enddo
             end if
          end do
       end do
    end do

    ! now, for each direction in turn, we need to apply a number of
    ! 'spreading operations'. Do z first... put blip(z) in 2, and
    ! del2blip(z) in 3

    call copy (blip_info(spec)%FullArraySize * this_nsf, work1, 1, work2, 1)
    call scal (blip_info(spec)%FullArraySize * this_nsf, FAC(0), work2, 1)
    do dz = 1, MAX_D
       offset = dz * blip_info(spec)%OneArraySize * &
            blip_info(spec)%OneArraySize * this_nsf
       call axpy ((blip_info(spec)%FullArraySize * this_nsf-offset), &
            FAC(dz), work1(1:), 1, work2(1+offset:), 1 )
       call axpy ((blip_info(spec)%FullArraySize * this_nsf-offset), &
            FAC(dz), work1(1+offset:), 1, work2(1:), 1 )
    end do

    ! now do y : put blip(y).blip(z) in 4,

    call copy (blip_info(spec)%FullArraySize * this_nsf, work2, 1, work4, 1)
    call scal (blip_info(spec)%FullArraySize * this_nsf, FAC(0), work4, 1)
    do dy = 1, MAX_D
       offset = dy * blip_info(spec)%OneArraySize * this_nsf
       call axpy ((blip_info(spec)%FullArraySize * this_nsf-offset), &
            FAC(dy), work2(1:), 1, work4(1+offset:), 1)
       call axpy ((blip_info(spec)%FullArraySize * this_nsf-offset), &
            FAC(dy), work2(1+offset:), 1, work4(1:), 1)
    end do

    ! and x - put it all into 6

    call copy (blip_info(spec)%FullArraySize * this_nsf, work4, 1, work6, 1)
    call scal (blip_info(spec)%FullArraySize * this_nsf, FAC(0), work6, 1)
    do dx = 1, MAX_D
       offset = dx * this_nsf
       call axpy((blip_info(spec)%FullArraySize * this_nsf-offset), &
            FAC(dx), work4(1:), 1, work6(1+offset:), 1 )
       call axpy((blip_info(spec)%FullArraySize * this_nsf-offset), &
            FAC(dx), work4(1+offset:), 1, work6(1:), 1 )
    end do

    ! now accumulate work6 onto the gradient
    call scal(blip_info(spec)%FullArraySize * this_nsf, &
              -two * spin_factor *                      &
              blip_info(spec)%SupportGridSpacing *      &
              blip_info(spec)%SupportGridSpacing *      &
              blip_info(spec)%SupportGridSpacing,       &
              work6, 1)

    offset = blip_info(spec)%BlipArraySize + 1
    do dx = - blip_info(spec)%BlipArraySize, blip_info(spec)%BlipArraySize
       do dy = - blip_info(spec)%BlipArraySize, blip_info(spec)%BlipArraySize
          do dz = - blip_info(spec)%BlipArraySize, blip_info(spec)%BlipArraySize
             l = blip_info(spec)%blip_number(dx,dy,dz)
             if (l .ne. 0) then
                at = (((dz + offset) * blip_info(spec)%OneArraySize + &
                     (dy + offset)) * blip_info(spec)%OneArraySize + &
                     (dx + offset)) * this_nsf
                do nsf1 = 1, this_nsf
                   this_blip_grad%supp_func(nsf1)%coefficients(l) = &
                        this_blip_grad%supp_func(nsf1)%coefficients(l) + &
                        work6(nsf1 + at)
                enddo
             end if
          end do
       end do
    end do
    call start_timer (tmr_std_allocation)
    deallocate (work1, work2, work4, work6, STAT=stat)
    if (stat /= 0) &
         call cq_abort ("Error deallocating arrays for onsite S blip grad: ", &
         blip_info(spec)%FullArraySize, this_nsf)
    call reg_dealloc_mem (area_minE, 6 * blip_info(spec)%FullArraySize * &
         this_nsf, type_dbl)
    call stop_timer (tmr_std_allocation)
    return
  end subroutine get_onsite_S_gradient
  !!***

end module blip_gradient
