! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module H_matrix_module
! ------------------------------------------------------------------------------
! Code area 3: Operators
! ------------------------------------------------------------------------------

!!****h* Conquest/H_matrix_module *
!!  NAME
!!   H_matrix_module
!!  USES
!!   atoms, blip_grid_transform_module, blip,
!!   calc_matrix_elements_module, common, cover_module, datatypes,
!!   dimens, generic_blas, generic_comms, grid_index, group_module,
!!   matrix_data, maxima_module, multiply_module, mult_module,
!!   numbers, primary_module, pseudopotential_data,
!!   set_blipgrid_module, set_bucket_module, species_module, workspace
!!  PURPOSE
!!   Collects together all routines needed for making the
!!   H matrix
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   06/06/2001
!!  MODIFICATION HISTORY
!!   11/06/2001 dave
!!    Added dependence on GenComms
!!   21/06/2001 dave
!!    Included get_onsite_T, get_xc_potential and matrix_scale
!!    Removed use matrix_scaling from get_HNL_matrix
!!   13/05/2002 dave
!!    Added headers, RCS object Id and a few more comments
!!   24/06/2002 dave
!!    Added comments and explicit NSF loop to matrix_scale
!!   29/07/2002 drb
!!    Changed get_HNL_matrix to use a local variable for transposes
!!   10:39, 04/02/2003 drb
!!    Small changes to get_HNL_matrix (use temp_TCS)
!!   11:51, 17/03/2003 drb
!!    Bug fix in get_T_matrix
!!   15:30, 08/04/2003 drb
!!    Added get_GTH_xc_potential as an alternative LDA parameterisation
!!   08:07, 2003/09/22 dave
!!    Added PAO/blip switches to get_HNL_matrix and get_T_matrix
!!   10:09, 13/02/2006 drb
!!    Removed all explicit references to data_ variables and rewrote
!!    in terms of new
!!    matrix routines
!!   2006/03/06 04:46 dave
!!    Removing passed arrays support, workspace_support
!!   2008/02/01 17:47 dave
!!    Changes for output to file not stdout
!!   2008/05/22 ast
!!    Added timers
!!   2012/03/15 L.Tong
!!    Changed spin implementation
!!  TODO
!!   08:28, 2003/09/22 dave
!!    Understand and document onsite_T
!   2014/09/15 18:30 lat
!!    fixed call start/stop_timer to timer_module (not timer_stdlocks_module !)
!   2019/12/30 tsuyoshi
!!    introduced flag_DumpChargeDensity to control dump_charge in the end of get_H_matrix
!!  2020/01/02 16:53 dave
!!    Moved flag to density_module.f90 (more logical location)
!!  SOURCE
!!
module H_matrix_module

  use global_module,          only: io_lun
  use timer_module,           only: start_timer, stop_timer, stop_print_timer
  use timer_module,           only: start_backtrace, stop_backtrace
  use timer_stdclocks_module, only: tmr_std_hmatrix,         &
                                    tmr_std_allocation,      &
                                    tmr_std_matrices

  implicit none

  ! Area identification
  integer, parameter, private :: area = 3

  logical :: locps_output
  integer :: locps_choice

!!***

contains

  ! ------------------------------------------------------------------------------
  ! Subroutine get_H_matrix
  ! ------------------------------------------------------------------------------

  !!****f* H_matrix_module/get_H_matrix_nospin *
  !!
  !! NAME
  !!  get_H_matrix
  !! USAGE
  !!
  !! PURPOSE
  !!  gets the entire Hamiltonian matrix
  !! INPUTS
  !!
  !!
  !! USES
  !!
  !! AUTHOR
  !!  C.M.Goringe
  !! CREATION DATE
  !!  11/10/95
  !! MODIFICATION HISTORY
  !!  10/08/2000 TM
  !!   New get_matrix_elements
  !!  14/05/2001 DRB
  !!   Added ROBODoc header and removed unnecessary arguments
  !!  17/05/2001 dave
  !!   Reduced call to get_T_matrix
  !!   Reduced call to get_h_on_support
  !!   Reduced call to get_HNL_matrix
  !!  18/05/2001 dave
  !!   Reduced overall argument list
  !!  07/06/2001 dave
  !!   Included in H_matrix_module
  !!  11/06/2001 dave
  !!   Changed to use GenComms for gsum
  !!  13/05/2002 dave
  !!   Added more comments and shifted format statements around
  !!  15/11/2002 tm
  !!   Changed pseudopotential_data use to pseudopotential_common
  !!  09:11, 2003/03/11 dave
  !!   Removed find_chdens argument and writing out of energies and
  !!   added flag to force or not the rebuilding of data_NL and data_KE
  !! (the non-local and kinetic energy contributions to the H matrix) 
  !!  2008/05/22 ast
  !!   Added timers
  !!  2011/10/06 13:51 dave
  !!   Added cDFT call to store H without constraint
  !!  2011/11/28 L.Tong
  !!   Removed redundant dependence on module global variable potential
  !!   from potential_module
  !!  2012/02/14 L.Tong
  !!   Renamed the subroutine to get_H_matrix_nospin, to be used with
  !!   the get_H_matrix interface
  !!  2012/03/15 L.Tong
  !!  - Major rewrite of the spin implementation
  !!  - Merged get_H_matrix_nospin and get_H_matrix_spin into one
  !!    subroutine
  !!  - Renamed the new subroutine get_H_matrix as it used to be.
  !!    Deleted het_H_matrix interface.
  !!  2012/04/16 L.Tong
  !!  - Moved all the XC functionals and related subroutines to XC_module
  !!  2012/04/26 16:13 dave
  !!   Changes for analytic evaluation of KE
  !!  2013/07/10 11:23 dave
  !!   Bug fix for sum over two components of rho even without spin (and moved rho_total alloc/dealloc)
  !!  2015/06/08 lat 
  !!   Added EXX+spin and experimental backtrace 
  !!  2016/07/13 18:30 nakata
  !!   Renamed subroutine get_h_on_support -> get_h_on_atomfns
  !!   Renamed H_on_supportfns -> H_on_atomfns
  !!  2016/07/15 18:30 nakata
  !!   Renamed sf_H_sf_rem -> atomf_H_atomf_rem
  !!  2016/08/01 17:30 nakata
  !!   Use atomf_H_atomf_rem instead of pao_H_sf_rem
  !!  2016/08/08 15:30 nakata
  !!   Renamed supportfns -> atomfns
  !!  2016/09/16 21:30 nakata
  !!   Introduced matKEatomf, matNLatomf, matXatomf, matHatomf
  !!   Added matHatomf -> matH transformation (AtomF_to_SF_transform)
  !!  2016/12/19 18:30 nakata
  !!   Removed unused fn_on_grid, length, paolength
  !!   Removed dHrange, matdH, paof and sf, which are no longer needed.
  !!  2017/01/26 20:00 nakata
  !!   Added optional passed variables build_AtomF_matrix and transform_AtomF_to_SF
  !!  2018/01/30 10:06 dave
  !!   Moved call to NA projector matrix build inside the rebuild_KE_NL loop to improve
  !!   efficiency (it should have been in there in the first place)
  !!  2018/11/13 17:30 nakata
  !!   Changed matS, matKE, matNL and matNA to be spin_SF dependent
  !!  2019/01/31 16:00 nakata
  !!   Moved dump_matrix(NSmatrix) to sub:get_S_matrix
  !!  2021/07/19 15:46 dave
  !!   Removed writing out of charge density
  !! SOURCE
  !!
  subroutine get_H_matrix(rebuild_KE_NL, fixed_potential, electrons, &
                          rho, size, level, build_AtomF_matrix, transform_AtomF_to_SF)

    use datatypes
    use numbers 
    use matrix_data,                 only: Hrange, Srange
    use mult_module,                 only: matNL, matKE, matH,          &
                                           matrix_scale, matrix_sum,    &
                                           matrix_product,              &
                                           matS, matX, matNAatomf, matNA, &
                                           matNLatomf, matKEatomf,      &
                                           matHatomf, matXatomf,        &
                                           AtomF_to_SF_transform, &
                                           allocate_temp_matrix, free_temp_matrix, &
                                           S_trans, matrix_product_trace
    use matrix_data, only: mat, aHa_range
    use pseudopotential_common,      only: non_local, pseudopotential, flag_neutral_atom_projector
    use set_bucket_module,           only: rem_bucket, atomf_H_atomf_rem
    use calc_matrix_elements_module, only: get_matrix_elements_new
    use GenComms,                    only: gsum, end_comms,             &
                                           my_barrier, inode, ionode,   &
                                           cq_abort
    use global_module,               only: iprint_ops, iprint_exx,      &
                                           flag_basis_set, PAOs,        &
                                           atomf, sf,                   &
                                           IPRINT_TIME_THRES1,          &
                                           iprint_SC,                   &
                                           flag_perform_cDFT,           &
                                           area_ops, nspin, nspin_SF,   &
                                           spin_factor, blips,          &
                                           flag_analytic_blip_int,      &
                                           flag_neutral_atom, min_layer
    use functions_on_grid,           only: atomfns, H_on_atomfns,       &
                                           gridfunctions
    use io_module,                   only: dump_matrix, dump_blips,     &
                                           write_matrix
    use dimens,                      only: n_my_grid_points
    use memory_module,               only: reg_alloc_mem,               &
                                           reg_dealloc_mem, type_dbl
    use group_module,                only: parts
    use primary_module,              only: bundle
    use cover_module,                only: BCS_parts
    use timer_module,                only: cq_timer, WITH_LEVEL
    ! This is used to declare a local timer
    use cdft_data,                   only: matHzero,                    &
                                           cDFT_NumberAtomGroups,       &
                                           cDFT_Vc, matWc
!****lat<$
    use global_module,               only: flag_exx, exx_niter, exx_siter, exx_alpha
    use exx_kernel_default,          only: get_X_matrix
    use exx_module,                  only: get_X_params
    use exx_types,                   only: exx_hgrid, exx_psolver, exx_radius
    use exx_io,                      only: exx_global_write
!****lat>$
    use energy, only: local_ps_energy
    use density_module, only: flag_DumpChargeDensity
    use io_module,      only: dump_charge
    
    implicit none

    ! Passed variables
    integer, optional :: level
    integer :: size
    logical :: rebuild_KE_NL, fixed_potential
    real(double), dimension(:)   :: electrons
    real(double), dimension(:,:) :: rho
    logical, optional :: build_AtomF_matrix
    logical, optional :: transform_AtomF_to_SF

    ! local variables
    real(double), dimension(:), allocatable :: rho_total
    real(double)   :: kinetic_energy, nl_energy
    integer        :: stat, spin, spin_SF
    type(cq_timer) :: tmr_l_hmatrix
    type(cq_timer) :: backtrace_timer
    integer        :: backtrace_level
    character(len=9), dimension(2) :: print_exxspin
    logical        :: flag_build_Hatomf, flag_do_SFtransform

    print_exxspin(1) = 'spin up  '
    print_exxspin(2) = 'spin down'

!****lat<$
    if (       present(level) ) backtrace_level = level+1
    if ( .not. present(level) ) backtrace_level = -10
    call start_backtrace(t=backtrace_timer,who='get_H_matrix',where=area,&
         level=backtrace_level,echo=.true.)
!****lat>$

    ! timer
    call start_timer(tmr_std_hmatrix)            ! Total
    call start_timer(tmr_l_hmatrix, WITH_LEVEL)  ! Just this call

    stat = 0
    if (inode == ionode .and. iprint_ops > 3) &
         write (io_lun, fmt='(10x,"Entering get_H_matrix")')

    flag_build_Hatomf = .true.
    if (flag_basis_set == PAOs .and. present(build_AtomF_matrix)) flag_build_Hatomf = build_AtomF_matrix

    if (flag_build_Hatomf) then
       ! zero the H matrix (in Conquest format)
       do spin = 1, nspin
          call matrix_scale(zero, matHatomf(spin))
          ! zero H_on_atomfns
          gridfunctions(H_on_atomfns(spin))%griddata = zero
       end do
       !
       !
       if (rebuild_KE_NL) then
          if (inode == ionode .and. iprint_ops > 3)&
               & write(io_lun, fmt='(2x,"Rebuilding KE")')
          ! both matKEatomf and matNLatomf are independent of spin (only XC is spin dependent)
          if(.NOT.flag_analytic_blip_int.OR.flag_basis_set/=blips) &
               call matrix_scale(zero, matKEatomf)
          call matrix_scale(zero, matNLatomf)
          ! get the T matrix and the kinetic energy...
          if(.NOT.flag_analytic_blip_int.OR.flag_basis_set/=blips) &
               call get_T_matrix
          ! now, we do the non-local part (if we are doing it)
          if (non_local) then
             if (inode == ionode .and. iprint_ops > 3) &
                  write (io_lun, fmt='(2x,"Rebuilding NL")')
             call get_HNL_matrix
          end if
          if (iprint_ops > 4) call dump_matrix("NNL_atomf", matNLatomf, inode)
          if (iprint_ops > 4) call dump_matrix("NKE_atomf", matKEatomf, inode)
          if(flag_neutral_atom_projector) then
             call matrix_scale(zero, matNAatomf)
             call get_HNA_matrix(matNAatomf)
             if (iprint_ops > 4) call dump_matrix("NNA_atomf", matNAatomf, inode)
          end if
       end if
       !
       !
       ! from here on, workspace support becomes h_on_atomfns...
       ! in fact, what we are getting here is (H_local - T) acting on support
       call get_h_on_atomfns(iprint_ops + min_layer, fixed_potential, electrons, rho, size)
       !
       !
       if (inode == ionode .and. iprint_ops > 3) &
            write (io_lun, *) 'Doing integration'
       ! Do the integration - support holds <phi| and workspace_support
       ! holds H|phi>. Inode starts from 1, and myid starts from 0. 
       ! get_matrix_elements_new takes myid
       do spin = 1, nspin
          call get_matrix_elements_new(inode-1, rem_bucket(atomf_H_atomf_rem), &
                                       matHatomf(spin), atomfns, &
                                       H_on_atomfns(spin))
       end do
       if (inode == ionode .and. iprint_ops > 3) write (io_lun, *) 'Done integration'
       !
       !
       if (iprint_ops > 4) then
          if (nspin == 1) then
             call dump_matrix("Nl_atomf",    matHatomf(1), inode)
          else
             call dump_matrix("Nl_up_atomf", matHatomf(1), inode)
             call dump_matrix("Nl_dn_atomf", matHatomf(2), inode)
          end if
       end if

       ! add the kinetic energy and non-local matrices to give the
       ! complete H matrix
       do spin = 1, nspin
          call matrix_sum(one, matHatomf(spin), half, matKEatomf)
          call matrix_sum(one, matHatomf(spin), one,  matNLatomf)
          if(flag_neutral_atom_projector) call matrix_sum(one, matHatomf(spin), one, matNAatomf)
       end do
       !
       !
!****lat<$
       if (flag_exx) then
          ! Ugly stuff but for now that's ok. Purpose is to adapt EXX accuracy to
          ! the SCF covergence: closer to convergence finest is the grid
          !
          !if (inode==ionode) print*, 'exx_pulay_r0 = ', exx_pulay_r0
          if  ( exx_niter < exx_siter ) then
             ! For first H building use pure DFT. To be improved for Hartree-Fock
             if (inode == ionode .and. iprint_exx > 3) &
                  write (io_lun, *) 'EXX: first guess from DFT'
             !
          else
             !
             if (inode == ionode .and. iprint_exx > 3) &
                  write (io_lun, *) 'EXX: setting get_X_matrix'
             call get_X_params(backtrace_level)

             call exx_global_write() 
             !
             !if (inode == ionode .and. iprint_exx > 3) &
                  !write (io_lun, *) 'EXX: doing get_X_matrix'
             do spin = 1, nspin
                if (inode == ionode .and. iprint_exx > 3) &
                write (io_lun, *) 'EXX: doing get_X_matrix: ', print_exxspin(spin)
                call get_X_matrix(spin,backtrace_level)
             end do
             !
             !if (inode == ionode .and. iprint_exx > 3) &
                  !write (io_lun, *) 'EXX: done get_X_matrix'
             do spin = 1, nspin
                if (inode == ionode .and. iprint_exx > 3) &
                write (io_lun, *) 'EXX: done get_X_matrix: ', print_exxspin(spin)
                call matrix_sum(one, matHatomf(spin),-exx_alpha*half, matXatomf(spin)) 
             end do
             !
          end if
          !
          exx_niter = exx_niter + 1
          !
       end if
!****lat>$
    endif ! flag_build_Hatomf

    ! For PAO-based contracted SFs, atomic functions are PAOs
    !     so we should transform H from atomic-function basis to SF basis.
    ! For blips and one_to_one PAOs, atomic functions are SFs
    !     so no transformation is needed.
    if (atomf.ne.sf) then
       flag_do_SFtransform = .true.
       if (present(transform_AtomF_to_SF)) flag_do_SFtransform = transform_AtomF_to_SF
    else
       flag_do_SFtransform = .false.
    endif

    if (flag_do_SFtransform) then
       do spin_SF = 1, nspin_SF
          call AtomF_to_SF_transform(matKE(spin_SF), matKEatomf, spin_SF, Hrange)   ! only to output kinetic energy
          call AtomF_to_SF_transform(matNL(spin_SF), matNLatomf, spin_SF, Hrange)   ! only to output non-local PP energy
          if(flag_neutral_atom_projector) &
               call AtomF_to_SF_transform(matNA(spin_SF), matNAatomf, spin_SF, Hrange)   ! only to output neutral atom energy
       enddo
       do spin = 1, nspin
          if (flag_exx) call AtomF_to_SF_transform(matX(spin), matXatomf(spin), spin, Hrange)   ! only to output EXX energy
          call AtomF_to_SF_transform(matH(spin), matHatomf(spin), spin, Hrange)   ! total electronic Hamiltonian
       enddo
    endif
    !
    !
    ! dump matrices if required
    if (iprint_ops > 3) then
       if (nspin == 1) then
          call dump_matrix("NH",    matH(1), inode)
       else
          call dump_matrix("NH_up", matH(1), inode)
          call dump_matrix("NH_dn", matH(2), inode)
       end if
    end if
    if (iprint_ops > 4) then
       if (nspin_SF == 1) then
          call dump_matrix("NNL", matNL(1), inode)
          call dump_matrix("NKE", matKE(1), inode)
       else
          call dump_matrix("NNL_up", matNL(1), inode)
          call dump_matrix("NNL_dn", matNL(2), inode)
          call dump_matrix("NKE_up", matKE(1), inode)
          call dump_matrix("NKE_dn", matKE(2), inode)
       end if
    endif
    !
    !
    ! dump charges if required
    if (flag_DumpChargeDensity .or. iprint_SC > 3) then
       if(nspin==1) then
          allocate(rho_total(size), STAT=stat)
          if (stat /= 0) call cq_abort("Error allocating rho_total: ", size)
          call reg_alloc_mem(area_ops, size, type_dbl)
          rho_total(:) = spin_factor * rho(:,1)
          call dump_charge(rho_total, size, inode, spin=0)
          deallocate(rho_total, STAT=stat)
          if (stat /= 0) call cq_abort("Error deallocating rho_total")
          call reg_dealloc_mem(area_ops, size, type_dbl)
       else if (nspin == 2) then
          call dump_charge(rho(:,1), size, inode, spin=1)
          call dump_charge(rho(:,2), size, inode, spin=2)
       end if
    end if
    !
    !
    ! Store the new H matrix for cDFT
    if (flag_perform_cDFT) then
       do spin = 1, nspin
          call matrix_sum(zero, matHzero(spin), one, matH(spin))
       end do
    endif
    !
    !
    ! timer
    call stop_print_timer(tmr_l_hmatrix, "get_H_matrix", IPRINT_TIME_THRES1)
    call stop_timer(tmr_std_hmatrix)


!****lat<$
    call stop_backtrace(t=backtrace_timer,who='get_H_matrix',echo=.true.)
!****lat>$

    return
  end subroutine get_H_matrix
  !!***


  ! -----------------------------------------------------------
  ! Subroutine get_h_on_atomfns
  ! -----------------------------------------------------------

  !!****f* H_matrix_module/get_h_on_atomfns_nospin *
  !!
  !!  NAME
  !!   get_h_on_atomfns
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Generates H|phi> on the grid
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   E.H.Hernandez
  !!  CREATION DATE
  !!   15/3/95
  !!  MODIFICATION HISTORY
  !!   18/01/2001 TM
  !!    Added new set_grid modules
  !!   15/05/2001 dave
  !!    Changed call to get_electronic_density, F90 mode and
  !!    added ROBODoc header
  !!   17/05/2001 dave
  !!    Simplified call to get_xc_potential
  !!   07/06/2001 dave
  !!    Included in H_matrix_module
  !!   11/06/2001 dave
  !!    Changed to use GenComms for gsum
  !!   09:17, 2003/03/11 dave
  !!    Removed find_chdens and total_energy arguments and call to
  !!    get_electronic_density
  !!   2004/10/05 drb
  !!    Added hartree_module use
  !!   13:58, 2006/07/17
  !!    Added functional type selector
  !!   2008/04/02  M. Todorovic
  !!    Added local potential output
  !!   2010/11/10  TM
  !!    Moved the parts of the calculation of delta_E_xc from
  !!    subroutines get_**_xc_potential to get_h_on_support for the
  !!    implementation of PCC (partial core correction).
  !!   2011/03/31  M.Arita
  !!    Added the contribution from P.C.C.
  !!   Saturday, 2011/07/30 L.Tong
  !!    Moved the dependence of H_on_supportfns from using module to
  !!    input parameters as hOnSupportFns.
  !!   2012/02/14 L.Tong
  !!    Renamed get_h_on_support_nospin and to be used with interface
  !!    get_h_on_support
  !!   2012/03/15 L.Tong
  !!   - Major rewrite for spin implementation. Now spin and non-spin
  !!     calculations are merged into one subroutine
  !!   - Renamed back to get_h_on_support
  !!   - removed input parameter hOnSupportFns, use module variable
  !!     H_on_supportfns from functions_on_grid instead.
  !!   2012/05/29 L.Tong
  !!   - Changed output level and format of electron number information
  !!   2012/05/29 L.Tong
  !!   - Cleaned up xc-functonal selector. Now spin and non-spin
  !!     calculations share the same calls more or less.
  !!   2013/07/10 11:23 dave
  !!     Bug fix for sum over two components of rho even without spin
  !!   2014/09/24 L.Truflandier
  !!   - Added temporary PBE0 and HF
  !!   - optional output of x_energy only
  !!   2015/11/24 08:38 dave
  !!    Adjusted name of hartree_energy to hartree_energy_total_rho for neutral atom implementation
  !!   2016/07/13 18:30 nakata
  !!    Renamed subroutine get_h_on_support -> get_h_on_atomfns
  !!    Renamed H_on_supportfns -> H_on_atomfns
  !!   2016/07/20 16:30 nakata
  !!    Renamed naba_atm -> naba_atoms_of_blocks
  !!   2016/08/01 17:30 nakata
  !!    Introduced atomf instead of sf and paof
  !!   2016/08/08 15:30 nakata
  !!    Renamed supportfns -> atomfns
  !!   2016/12/19 18:30 nakata
  !!    Removed unused fn_on_grid
  !!   2018/02/13 12:17 dave
  !!    New XC interface implemented
  !!  SOURCE
  !!
  subroutine get_h_on_atomfns(output_level, fixed_potential, &
                              electrons, rho, size)

    use datatypes
    use numbers
    use global_module,               only: atomf,&
                                           nspin, spin_factor,         &
                                           flag_pcc_global, area_ops,  &
                                           exx_alpha, exx_niter, exx_siter, &
                                           flag_neutral_atom
     
    use XC,                          only: get_xc_potential
    use GenBlas,                     only: copy, axpy, dot, rsum
    use dimens,                      only: grid_point_volume,          &
                                           n_my_grid_points, n_grid_z

    use block_module,                only: n_blocks, n_pts_in_block
    use primary_module,              only: domain
    use set_blipgrid_module,         only: naba_atoms_of_blocks
    use density_module,              only: density_pcc, density_atom
    use GenComms,                    only: gsum, inode, ionode, cq_abort
    use energy,                      only: hartree_energy_total_rho,  &
                                           xc_energy,       &
                                           x_energy,        &
                                           local_ps_energy, &
                                           delta_E_hartree, &
                                           hartree_energy_drho, &
                                           hartree_energy_drho_atom_rho, &
                                           delta_E_xc
    use hartree_module,              only: hartree, hartree_stress
    use functions_on_grid,           only: gridfunctions, &
                                           atomfns, H_on_atomfns

    use calc_matrix_elements_module, only: norb
    use pseudopotential_common,      only: pseudopotential, flag_neutral_atom_projector
    use potential_module,            only: potential
    use maxima_module,               only: maxngrid
    use memory_module,               only: reg_alloc_mem,              &
                                           reg_dealloc_mem, type_dbl
    use fft_module,                  only: fft3, hartree_factor,       &
                                           z_columns_node, i0
    use io_module,                   only: dump_locps

    implicit none

    ! Passed variables
    integer, intent(in) :: output_level, size
    logical, intent(in) :: fixed_potential
    real(double), dimension(:,:), intent(in)  :: rho
    real(double), dimension(:),   intent(out) :: electrons

    ! Local variables
    integer :: n, m, nb, atom, nsf1, point, stat, i, pot_flag, igrid, spin, exx_nit
    real(double) :: fften, electrons_tot, exx_tmp, temp_stress
    logical     , dimension(4)    :: dump_pot
    real(double), dimension(:),   allocatable :: xc_epsilon ! energy_density of XC
    real(double), dimension(:),   allocatable :: h_potential
    real(double), dimension(:),   allocatable :: rho_tot
    real(double), dimension(:,:), allocatable :: xc_potential
    real(double), dimension(:,:), allocatable :: density_wk ! rho + density_pcc
    real(double), dimension(:),   allocatable :: density_wk_tot
    !complex(double_cplx), dimension(:), allocatable :: chdenr, locpotr
    real(double), dimension(:),   allocatable :: drho_tot

    ! for Neutral atom potential
    if( flag_neutral_atom ) then
       allocate( drho_tot(size), STAT=stat)
       if (stat /= 0) &
            call cq_abort("get_h_on_atomfns: Error allocating drho_tot:", size )
       call reg_alloc_mem(area_ops, size, type_dbl)
    end if

    
    allocate(xc_epsilon(size), h_potential(size), rho_tot(size), &
             xc_potential(size,nspin), STAT=stat)
    if (stat /= 0) &
         call cq_abort("get_h_on_atomfns: Error allocating mem:", size, nspin)
    call reg_alloc_mem(area_ops, (3+nspin)*size, type_dbl)
    
    !call dump_locps(pseudopotential,size,inode)
    !allocate(chdenr(size), locpotr(size), STAT=stat)
    !call fft3( rho, chdenr, size, -1 )
    !call fft3( pseudopotential, locpotr, size, -1 )
    !fften = zero
    !do i = 1, z_columns_node(inode)*n_grid_z
    !   if(i==i0) then
    !      write (io_lun, *) 'G=0: ', &
    !                        grid_point_volume * real(chdenr(i),double) * &
    !                        real(locpotr(i),double), &
    !                        aimag(chdenr(i)) * aimag(locpotr(i)) * &
    !                         grid_point_volume
    !      !locpotr(i0)=cmplx(zero,zero)
    !   else
    !      fften = fften + real(chdenr(i),double) * real(locpotr(i),double) * &
    !              grid_point_volume - aimag(chdenr(i)) * aimag(locpotr(i)) * &
    !              grid_point_volume
    !   end if
    !end do
    !write (io_lun,*) 'Energy via FFT: ', fften
    !call fft3(pseudopotential, locpotr, size, 1)
    !deallocate(chdenr, locpotr, STAT=stat)
    !
    !
    ! first initialise some arrays
    h_potential = zero
    do spin = 1, nspin
       gridfunctions(H_on_atomfns(spin))%griddata = zero
       potential(:,spin)    = zero
       xc_potential(:,spin) = zero
    end do
    rho_tot = zero
    do spin = 1, nspin
       rho_tot(:) = rho_tot(:) + spin_factor*rho(:,spin)
       electrons(spin) = grid_point_volume * rsum(n_my_grid_points, rho(:,spin), 1)
       call gsum(electrons(spin))
    end do
    ! for Neutral atom potential
    if( flag_neutral_atom ) then
       drho_tot(:) = rho_tot(:) - density_atom(:)
    end if
    
    electrons_tot = spin_factor * sum(electrons(:))
    if (inode == ionode .and. output_level >= 3) then
       write (io_lun, '(10x,a)') &
            'get_h_on_atomfns: Electron Count, up, down and total:'
       write (io_lun, '(10x, 3f25.15)') &
            electrons(1), electrons(nspin), electrons_tot
    end if
    !
    !
    ! now calculate the hartree potential on the grid
    if(flag_neutral_atom) then
       call hartree(drho_tot, h_potential, maxngrid, hartree_energy_drho)
       ! At this point, hartree_stress should be correct
       ! get the pseudopotential energy
       hartree_energy_drho_atom_rho = &
            dot(n_my_grid_points, h_potential, 1, density_atom, 1) * &
            grid_point_volume
       call gsum(hartree_energy_drho_atom_rho)
       ! This is the effective correction term
       delta_E_hartree = - hartree_energy_drho - hartree_energy_drho_atom_rho
    else
       call hartree(rho_tot, h_potential, maxngrid, hartree_energy_total_rho)
       !
       !
       ! correction term
       delta_E_hartree = - hartree_energy_total_rho
    end if
    !
    !
    ! for P.C.C.
    if (flag_pcc_global) then
       allocate(density_wk(size,nspin), density_wk_tot(size), STAT=stat)
       if (stat /= 0) &
            call cq_abort("Error allocating density_wk, density_wk_tot: ", &
                          stat)
       call reg_alloc_mem(area_ops, size * (nspin + 1), type_dbl)
       !density_wk = zero
       do spin = 1, nspin
          density_wk(:,spin) = rho(:,spin) + half * density_pcc(:)
       end do
       density_wk_tot = rho_tot + density_pcc
    end if
    !  
    !
    if (flag_pcc_global) then
       call get_xc_potential(density=density_wk, size=size, &
            xc_potential=xc_potential,    &
            xc_epsilon  =xc_epsilon, &
            xc_energy   =xc_energy,  &
            x_energy    =x_energy    )
    else
       call get_xc_potential(density=rho, size=size,     &
            xc_potential=xc_potential, &
            xc_epsilon  =xc_epsilon,        & 
            xc_energy   =xc_energy,         &
            x_energy    =x_energy)
    end if
    !
    !
    ! Calculation of delta_E_xc
    delta_E_xc = zero
    if (flag_pcc_global) then
       do igrid = 1, n_my_grid_points
          delta_E_xc = delta_E_xc + &
                       xc_epsilon(igrid) * density_wk_tot(igrid) - &
                       xc_potential(igrid,1) * rho(igrid,1) - &
                       xc_potential(igrid,nspin) * rho(igrid,nspin)
       end do
    else
       do igrid = 1, n_my_grid_points
          delta_E_xc = delta_E_xc + &
                       xc_epsilon(igrid) * rho_tot(igrid) - &
                       xc_potential(igrid,1) * rho(igrid,1) - &
                       xc_potential(igrid,nspin) * rho(igrid,nspin)
       end do
    end if ! (flag_pcc_global)
    call gsum(delta_E_xc)
    delta_E_xc = delta_E_xc * grid_point_volume
    !
    !
    ! Make total potential
    if (.not. fixed_potential) then
       do spin = 1, nspin
          call copy(n_my_grid_points, h_potential, 1, potential(:,spin), 1)
          call axpy(n_my_grid_points, one, xc_potential(:,spin), 1, &
                    potential(:,spin), 1)
          if(.NOT.flag_neutral_atom_projector) &
               call axpy(n_my_grid_points, one, pseudopotential, 1, &
               potential(:,spin), 1)
       end do
    end if
    !
    !
    ! Print potential, if necessary
    if (locps_output) then
       dump_pot = (/.false., .true., .false., .false./)
       if (locps_choice < 0 .or. locps_choice > 15) then
          if (inode == ionode) &
               write (io_lun, *) 'Bad choice for local potential &
                                  &printout: no output.'
       else
          if (inode == ionode) &
               write (io_lun, *) 'Writing local potential to file(s).'
          do i = 4, 1, -1
             pot_flag = mod(locps_choice/(2**(i-1)), 2)
             if (pot_flag > 0) dump_pot(i) = .true.
             locps_choice = locps_choice - pot_flag*(2**(i-1))
          end do
       end if
       if (dump_pot(1)) call dump_locps("Hartree", h_potential, size, inode)
       if (dump_pot(2)) then
          if (nspin == 1) then
             call dump_locps("XC", xc_potential(:,1), size, inode)
          else
             call dump_locps("XC_up", xc_potential(:,1), size, inode)
             call dump_locps("XC_dn", xc_potential(:,2), size, inode)
          end if
       end if
       if (dump_pot(3)) call dump_locps("PS", pseudopotential, size, inode)
       if (dump_pot(4)) then
          if (nspin == 1) then
             call dump_locps("Total", potential(:,1), size, inode)
          else
             call dump_locps("Total_up", potential(:,1), size, inode)
             call dump_locps("Total_dn", potential(:,2), size, inode)
          end if
       end if
    end if
    !
    !
    ! get the pseudopotential energy
    ! In the OpenMX VNA formulation, we need the full integral over VNA and rhoSCF
    ! NB if we use blips or projectors for NA this should become dot(K,HNA)
    if(.NOT.flag_neutral_atom_projector) then
       local_ps_energy = &
            dot(n_my_grid_points, pseudopotential, 1, rho_tot, 1) * &
            grid_point_volume
       call gsum(local_ps_energy)
    end if
    !
    !
    ! now act with the potential on the support functions to get
    ! (H_local - T) on support
    n = 0
    m = 0
    do nb = 1, domain%groups_on_node
       if (naba_atoms_of_blocks(atomf)%no_of_atom(nb) > 0) then
          do atom = 1, naba_atoms_of_blocks(atomf)%no_of_atom(nb)
             do nsf1 = 1, norb(naba_atoms_of_blocks(atomf),atom,nb)
                do point = 1, n_pts_in_block
                   n = n + 1
                   do spin = 1, nspin
                      gridfunctions(H_on_atomfns(spin))%griddata(n) = &
                           gridfunctions(atomfns)%griddata(n) *  &
                           potential(m+point,spin)
                   end do ! spin
                end do ! point
             end do ! nsf1
          end do ! atom
       end if ! (naba_atoms_of_blocks(atomf)%no_of_atom(nb) > 0)
       m = m + n_pts_in_block
    end do ! nb
    !
    !
    if (flag_pcc_global) then
       deallocate(density_wk, STAT=stat)
       if (stat /= 0) &
            call cq_abort("Error deallocating density_wk: ", stat)
       call reg_dealloc_mem(area_ops, size, type_dbl)
    endif
    !
    !
    deallocate(xc_epsilon, h_potential, rho_tot, xc_potential, STAT=stat)
    if (stat /= 0) call cq_abort("get_h_on_atomfns: Error deallocating mem")
    call reg_dealloc_mem(area_ops, (3+nspin)*size, type_dbl)

    ! for Neutral atom potential
    if( flag_neutral_atom ) then
       deallocate(drho_tot, STAT=stat)
       if (stat /= 0) call cq_abort("get_h_on_atomfns: Error deallocating mem")
       call reg_dealloc_mem(area_ops, size, type_dbl)
    end if
    return
  end subroutine get_h_on_atomfns
  !!***

  ! -----------------------------------------------------------
  ! Subroutine get_HNL_matrix
  ! -----------------------------------------------------------

  !!****f* H_matrix_module/get_HNL_matrix *
  !!
  !!  NAME
  !!   get_HNL_matrix
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Generates the non-local pseudopotential part of the Hamiltonian matrix
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   25/4/96 CMG
  !!  MODIFICATION HISTORY
  !!   10/08/2000 TM
  !!    New get_matrix_elements
  !!   14/05/01 DRB
  !!    Added ROBODoc header and shortened passed variable list
  !!   17/05/2001 dave
  !!    Simplified passed variables and modules used
  !!   07/06/2001 dave
  !!    Included in H_matrix_module
  !!   11/06/2001 dave
  !!    Changed to use GenComms for gsum
  !!   21/06/2001 dave
  !!    Removed use matrix_scaling (now in this module)
  !!   29/07/2002 drb
  !!    Added a local transpose variable to correct a bug when NCF/=NSF
  !!   15/11/2002 tm
  !!    Changed references for pseudopotentials to pseudopotential_common
  !!    Added matrix_scale_tm as internal subroutine
  !!   10:38, 04/02/2003 drb
  !!    The transpose now uses temp_TCS (a global variable)
  !!   08:03, 2003/09/22 dave
  !!    Added blip/PAO switch; removed MAX_H_ELEMENTS and MAX_SC_ELEMENTS
  !!   08:29, 2004/07/23 dave
  !!    Added lines to allow calculation of derivative of
  !!    <\phi|\hat{NL}|\phi> with PAOs as basis for \phi
  !!   2011/08/25 L.Tong
  !!    Added Spin polarisation
  !!     note that the non-local potential is it-self independent of
  !!     spin and hence does not require modification. However matdH is
  !!     spin dependent due to matH is spin dependent. And hence the
  !!     contribution to matdH from the non local potential must also be
  !!     added to matdH_dn. The value of the contribution however is
  !!     same for both spin components
  !!   2016/07/15 18:30 nakata
  !!    Renamed sf_nlpf_rem -> atomf_nlpf_rem
  !!   2016/07/29 18:30 nakata
  !!    Renamed supports_on_atom -> blips_on_atom
  !!   2016/08/01 17:30 nakata
  !!    Introduced atomf instead of sf
  !!   2016/08/08 15:30 nakata
  !!    Renamed supportfns -> atomfns
  !!   2016/09/20 18:30 nakata
  !!    Introduced matAP, matPA, matNLatomf, APrange, AP_trans, AP_PA_aHa
  !!    instead of matSC, matCS, matNL     , SPrange, SP_trans, SP_PS_H
  !!    matNLatomf is used from module, not from passed variables
  !!   2016/12/19 18:30 nakata
  !!    Removed PAOPrange, dHrange, matdAP, matdCNL and matdH, 
  !!    which are no longer needed.
  !!  SOURCE
  !!
  subroutine get_HNL_matrix

    use datatypes
    use numbers
    use matrix_data, only: mat, APrange, halo
    use mult_module, only: mult, AP_PA_aHa, AP_trans,                  &
                           matPA, matAP, matNLatomf,                   &
                           allocate_temp_matrix,                       &
                           free_temp_matrix, matrix_product,           &
                           matrix_sum, matrix_transpose, matrix_scale
    use pseudopotential_data,        only: n_projectors, l_core, recip_scale
    use pseudopotential_common,      only: pseudo_type, OLDPS, SIESTA, &
                                           STATE, ABINIT
    use species_module,              only: species
    use set_bucket_module,           only: rem_bucket, atomf_nlpf_rem
    use calc_matrix_elements_module, only: get_matrix_elements_new
    use global_module,               only: flag_basis_set, PAOs,       &
                                           blips, nlpf, atomf,         &
                                           iprint_ops,                 &
                                           nspin, id_glob,             &
                                           species_glob,               &
                                           flag_analytic_blip_int, min_layer
    use GenComms,                    only: cq_abort, myid, inode, ionode
    use GenBlas,                     only: axpy
    use build_PAO_matrices,          only: assemble_2
    use functions_on_grid,           only: atomfns, pseudofns
    use io_module,                   only: dump_matrix
    use nlpf2blip,                   only: get_SP, nlpf_on_atom
    use primary_module ,             only: bundle
    use support_spec_format,         only: blips_on_atom
    use group_module,                only: parts
    use primary_module,              only: bundle
    use cover_module,                only: BCS_parts
    use species_module,              only: nsf_species, nlpf_species

    implicit none

    ! Passed variables
    integer :: matNL_tmp

    ! Local variables
    integer      :: stat, matAPtmp, np, ni, iprim, spec,                    &
                    this_nsf, this_nlpf, nab, ist, gcspart, n1, n2,         &
                    neigh_global_part, neigh_global_num, neigh_species, i1, &
                    i2, wheremat, spin
    real(double) :: dx, dy, dz


    matAPtmp = allocate_temp_matrix(APrange, AP_trans, atomf, nlpf)

    ! first, get the overlap of support functions with core
    ! pseudowavefunctions
    if (flag_basis_set == blips) then
       if (inode == ionode .and. iprint_ops > 2) &
            write(io_lun,*) 'Calling get_matrix_elements'
       if (flag_analytic_blip_int) then
          call matrix_scale(zero, matAP)
          iprim = 0
          do np = 1, bundle%groups_on_node
             if (bundle%nm_nodgroup(np) > 0) then
                do ni = 1, bundle%nm_nodgroup(np)
                   iprim = iprim + 1
                   spec = bundle%species(iprim)
                   ! write (60, *) "#Atom ", iprim
                   ! Loop over neighbours of atom
                   do nab = 1, mat(np,APrange)%n_nab(ni)
                      ist = mat(np,APrange)%i_acc(ni) + nab - 1
                      ! Build the distances between atoms - needed for phases
                      gcspart = &
                           BCS_parts%icover_ibeg(mat(np,APrange)%i_part(ist)) + &
                           mat(np,APrange)%i_seq(ist) - 1
                      ! Displacement vector
                      dx = BCS_parts%xcover(gcspart) - bundle%xprim(iprim)
                      dy = BCS_parts%ycover(gcspart) - bundle%yprim(iprim)
                      dz = BCS_parts%zcover(gcspart) - bundle%zprim(iprim)
                      ! We need to know the species of neighbour
                      neigh_global_part = &
                           BCS_parts%lab_cell(mat(np,APrange)%i_part(ist))
                      neigh_global_num  = &
                           id_glob(parts%icell_beg(neigh_global_part) + &
                                   mat(np,APrange)%i_seq(ist) - 1)
                      ! write (60, *) "#Nab no and glob: ", &
                      !               nab, neigh_global_num, dx, dy, dz
                      neigh_species = species_glob(neigh_global_num)
                      !write(io_lun,fmt='(2x,"Offset: ",3f7.2,4i4)') dx, dy, dz, iprim,neigh_global_num, spec,neigh_species
                      !do n1=1,nsf_species(spec)
                      !   do n2=1,nlpf_species(neigh_species)
                      !      call scale_matrix_value(matAP,np,ni,iprim,nab,n1,n2,zero)
                      !   end do
                      !end do
                      call get_SP(blips_on_atom(iprim),                  &
                                  nlpf_on_atom(neigh_species), matAP,    &
                                  iprim, halo(APrange)%i_halo(gcspart), &
                                  dx, dy, dz, spec, neigh_species)
                   end do
                   ! write (60, *) "&"
                end do
             end if
          end do
          !call dump_matrix("NSC2",matAP,inode)
       else
          call get_matrix_elements_new(myid, rem_bucket(atomf_nlpf_rem), &
                                       matAP, atomfns, pseudofns)
       end if
       !call dump_matrix("NSC",matAP,inode)
    else if (flag_basis_set == PAOs) then
       ! Use assemble to generate matrix elements
       if (inode == ionode .and. iprint_ops + min_layer > 3) &
            write (io_lun, *) 'Calling assemble'
       call assemble_2(APrange, matAP, 3)
       if (inode == ionode .AND. iprint_ops + min_layer > 3) &
            write(io_lun,*) 'Called assemble'
    else
       call cq_abort('get_HNL_matrix: basis set incorrectly specified ', &
                     flag_basis_set)
    end if
    if (inode == ionode .and. iprint_ops + min_layer > 3) write(io_lun,*) 'Made SP'

    call matrix_sum(zero, matAPtmp, one ,matAP)
    if (mult(AP_PA_aHa)%mult_type == 2) then ! type 2 means no transpose necessary
       select case (pseudo_type)
       case (OLDPS)
          call matrix_scale_diag(matAP, species, n_projectors, l_core,&
                                 recip_scale, APrange)
       case(SIESTA)
          call matrix_scale_diag_tm(matAP, APrange)
       case(ABINIT)
          call matrix_scale_diag_tm(matAP, APrange)
       end select
       call matrix_product(matAPtmp, matAP, matNLatomf, mult(AP_PA_aHa))
    else ! Transpose SP->PS, then mult
       if (inode == ionode .and. iprint_ops + min_layer > 3) &
            write (io_lun, *) 'Type 1 ', matAP, matPA
       call matrix_transpose(matAP, matPA)
       if (inode == ionode .and. iprint_ops + min_layer > 3) &
            write (io_lun, *) 'Done transpose'
       select case (pseudo_type)
       case (OLDPS)
          call matrix_scale_diag(matAP, species, n_projectors, l_core,&
                                 recip_scale, APrange)
       case (SIESTA)
          if (inode == ionode .and. iprint_ops + min_layer > 3) &
               write (io_lun, *) 'Doing scale'
          call matrix_scale_diag_tm(matAP, APrange)
       case(ABINIT)
          if (inode == ionode .and. iprint_ops + min_layer > 3) &
               write (io_lun, *) 'Doing scale'
          call matrix_scale_diag_tm(matAP, APrange)
       end select
       call matrix_product(matAP, matPA, matNLatomf, mult(AP_PA_aHa))
    endif

    call free_temp_matrix(matAPtmp)
    return

  contains

    ! Added nonef to allow scaling of the matrix <PAO|projector>
    ! required for PAO basis of support functions
    ! DRB 2004/07/23 08:28
    subroutine matrix_scale_diag_tm(matAP, range)

      ! Module usage
      use datatypes,      only: double
      use numbers,        only: RD_ERR, zero
      use group_module,   only: parts
      use primary_module, only: bundle
      use cover_module,   only: BCS_parts
      use matrix_data,    only: mat
      use species_module, only: species
      use pseudo_tm_info, only: pseudo
      use mult_module,    only: scale_matrix_value, store_matrix_value
      use GenComms,       only: inode, ionode

      implicit none

      ! Passed variables
      integer :: range, matAP

      ! Local variables
      integer :: np, i, nb, isu, ind_cover, ind_qart, ip
      integer :: k, species_k, nl, nsf1
      integer :: mmax, mm, n_proj

      ! Loop over the i elements and apply the appropriate scaling factor
      ip = 0
      call start_timer(tmr_std_matrices)
      do np = 1,bundle%groups_on_node
         if (bundle%nm_nodgroup(np) > 0) then
            do i = 1,bundle%nm_nodgroup(np)
               ip = ip+1
               if (mat(np,range)%n_nab(i) > 0) then
                  do nb = 1, mat(np,range)%n_nab(i)
                     isu = mat(np,range)%i_acc(i)+nb-1
                     ind_cover=mat(np,range)%i_part(isu)
                     ind_qart=BCS_parts%lab_cell(ind_cover)
                     k = parts%icell_beg(ind_qart)+mat(np,range)%i_seq(isu)-1
                     species_k = species(k)
                     n_proj=0
                    if (pseudo(species_k)%n_pjnl > 0) then
                     do nl = 1, pseudo(species_k)%n_pjnl
                        mmax=pseudo(species_k)%pjnl_l(nl)*2+1
                        do mm=1,mmax
                           n_proj=n_proj+1
                           do nsf1=1,mat(np,range)%ndimi(i)
                              call scale_matrix_value(&
                                   matAP, np, i, ip, nb, nsf1, n_proj, &
                                   pseudo(species_k)%pjnl_ekb(nl))
                              if (abs(pseudo(species_k)%pjnl_ekb(nl)) < &
                                  RD_ERR .and. inode == ionode) &
                                   write (io_lun, *) &
                                   'ekb = 0!!   for nl, species_k, ekb = ', &
                                   nl, species_k, pseudo(species_k)%pjnl_ekb(nl)
                           enddo !nsf1=1,nonef
                        enddo !mm=1,mmax
                     enddo !nl = 1, pseudo(species_k)%n_pjnl
                    endif ! (pseudo(species_k)%n_pjnl > 0
                  enddo ! nb=mat%n_nab(i)
               endif ! mat%n_nab(i)>0
            enddo ! i=bundle%nm_nodgroup
         endif ! bundle%nm_nodgroup>0
      enddo ! np=bundle%groups_on_node
      call stop_timer(tmr_std_matrices)
      return
    end subroutine matrix_scale_diag_tm

  end subroutine get_HNL_matrix
  !!***

  ! -----------------------------------------------------------
  ! Subroutine get_HNA_matrix
  ! -----------------------------------------------------------

  !!****f* H_matrix_module/get_HNA_matrix *
  !!
  !!  NAME
  !!   get_HNA_matrix
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Generates the neutral-atom matrix using PAO integrals for 
  !!   1- and 2-centre terms, and projectors for the 3-centre terms
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   2017/11/10
  !!  MODIFICATION HISTORY
  !!   2018/01/25 12:50 JST dave
  !!    Changed transpose type for matNA
  !!  SOURCE
  !!
  subroutine get_HNA_matrix(matNA)

    use datatypes
    use numbers
    use matrix_data, only: mat, aNArange, aHa_range, halo
    use mult_module, only: matrix_sum, matrix_transpose, store_matrix_value, return_matrix_value, &
         allocate_temp_matrix, free_temp_matrix, S_trans, scale_matrix_value,  aNA_NAa_aHa, &
         aNA_trans, matNAa, mataNA, matrix_scale, mult, matrix_product, matrix_pos, return_matrix_len, matK, &
         matrix_product_trace, aNAa_trans
    use species_module,              only: species
    use global_module,               only: flag_basis_set, PAOs,       &
                                           blips, nlpf, atomf, iprint_ops, &
                                           nspin, id_glob,             &
                                           species_glob,               &
                                           flag_analytic_blip_int
    use build_PAO_matrices,          only: assemble_2
    use primary_module ,             only: bundle
    use pao_format, ONLY: pao
    use pseudo_tm_info, only: pseudo
    use cover_module,   only: BCS_parts
    use group_module,   only: parts
    use io_module, only: dump_matrix
    use GenComms, only: inode

    implicit none

    ! Passed variables
    integer :: matNA

    ! Local variables
    integer      :: matNAT, np, i, ip, nsf1, nsf2, ist, nb, gcspart, wheremat
    real(double) :: val

    ! 1- and 2-centre terms: <ii|j> and <i|jj> by assemble, transpose and sum
    matNAT = allocate_temp_matrix(aHa_range, aNAa_trans, atomf, atomf)
    call assemble_2(aHa_range, matNA, 4)
    call matrix_transpose(matNA, matNAT)
    call matrix_sum(one,matNA,one,matNAT)
    ! Scale on-site by half to avoid double-counting
    ip = 0
    !call start_timer(tmr_std_matrices)
    do np = 1,bundle%groups_on_node
       if (bundle%nm_nodgroup(np) > 0) then
          do i = 1,bundle%nm_nodgroup(np)
             ip = ip+1
             do nsf1=1,mat(np,aHa_range)%ndimi(i)
                do nsf2=1,mat(np,aHa_range)%ndimi(i)
                   call scale_matrix_value(matNA, np, i, ip, 0, nsf1, nsf2, half,1)
                end do
             end do
          end do
       end if
    end do
    ! < i | Vj | i > terms
    call matrix_scale(zero,matNAT)
    ! Assemble 2 accumulates on-site
    call assemble_2(aHa_range, matNAT, 6) 
    call matrix_sum(one,matNA,one,matNAT)
    ! 3-centre terms
    call assemble_2(aNArange, mataNA, 5)
    ! Transpose
    call matrix_transpose(mataNA, matNAa)
    ! Scale
    call matrix_scale_diag_NA(mataNA, aNArange)
    ! And form the matrix
    call matrix_scale(zero,matNAT)
    call matrix_product(mataNA, matNAa, matNAT, mult(aNA_NAa_aHa))
    ! Zero on-site terms
    ip=0
    do np = 1,bundle%groups_on_node
       if (bundle%nm_nodgroup(np) > 0) then
          do i = 1,bundle%nm_nodgroup(np)
             ip = ip+1
             do nsf1=1,mat(np,aHa_range)%ndimi(i)
                do nsf2=1,mat(np,aHa_range)%ndimi(i)
                   call scale_matrix_value(matNAT, np, i, ip, 0, nsf1, nsf2, zero,1)
                end do
             end do
          end do
       end if
    end do
    call matrix_sum(one,matNA,one,matNAT)
    call free_temp_matrix(matNAT)
    return

  contains
    subroutine matrix_scale_diag_NA(mataNA, range)

      ! Module usage
      use datatypes,      only: double
      use numbers,        only: RD_ERR, zero
      use group_module,   only: parts
      use primary_module, only: bundle
      use cover_module,   only: BCS_parts
      use matrix_data,    only: mat
      use species_module, only: species
      use pseudo_tm_info, only: pseudo
      use mult_module,    only: scale_matrix_value, return_matrix_value, store_matrix_value
      use GenComms,       only: inode, ionode

      implicit none

      ! Passed variables
      integer :: range, mataNA

      ! Local variables
      integer :: np, i, nb, isu, ind_cover, ind_qart, ip
      integer :: k, species_k, nl, nsf1
      integer :: mmax, mm, n_proj

      ! Loop over the i elements and apply the appropriate scaling factor
      ip = 0
      call start_timer(tmr_std_matrices)
      do np = 1,bundle%groups_on_node
         if (bundle%nm_nodgroup(np) > 0) then
            do i = 1,bundle%nm_nodgroup(np)
               ip = ip+1
               if (mat(np,range)%n_nab(i) > 0) then
                  do nb = 1, mat(np,range)%n_nab(i)
                     isu = mat(np,range)%i_acc(i)+nb-1
                     ind_cover=mat(np,range)%i_part(isu)
                     ind_qart=BCS_parts%lab_cell(ind_cover)
                     k = parts%icell_beg(ind_qart)+mat(np,range)%i_seq(isu)-1
                     species_k = species(k)
                     n_proj=0
                    if (pseudo(species_k)%n_pjna > 0) then
                     do nl = 1, pseudo(species_k)%n_pjna
                        mmax=pseudo(species_k)%pjna_l(nl)*2+1
                        do mm=1,mmax
                           n_proj=n_proj+1
                           do nsf1=1,mat(np,range)%ndimi(i)
                              call scale_matrix_value(&
                                   mataNA, np, i, ip, nb, nsf1, n_proj, &
                                   pseudo(species_k)%pjna_ekb(nl))
                           enddo !nsf1=1,nonef
                        enddo !mm=1,mmax
                     enddo !nl = 1, pseudo(species_k)%n_pjna
                  endif ! (pseudo(species_k)%n_pjna > 0
                  enddo ! nb=mat%n_nab(i)
               endif ! mat%n_nab(i)>0
            enddo ! i=bundle%nm_nodgroup
         endif ! bundle%nm_nodgroup>0
      enddo ! np=bundle%groups_on_node
      call stop_timer(tmr_std_matrices)
      return
    end subroutine matrix_scale_diag_NA
    
  end subroutine get_HNA_matrix
  
  ! -----------------------------------------------------------
  ! Subroutine get_T_matrix
  ! -----------------------------------------------------------

  !!****f* H_matrix_module/get_T_matrix *
  !!
  !!  NAME
  !!   get_T_matrix
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Calculates the kinetic energy matrix elements
  !!
  !!   The kinetic energy matrix is evaluated in two parts. For the blocks
  !!   off the leading diagonal (intersite terms) we integrate on the grid
  !!   Grad(psi_i).Grad(psi_j); this is done by looping over the three
  !!   cartesian components of the Grad and using blip_to_grad and
  !!   get_matrix_elements.
  !!
  !!  NAME
  !!   get_T_matrix
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Calculates the kinetic energy matrix elements
  !!
  !!   The kinetic energy matrix is evaluated in two parts. For the blocks
  !!   off the leading diagonal (intersite terms) we integrate on the grid
  !!   Grad(psi_i).Grad(psi_j); this is done by looping over the three
  !!   cartesian components of the Grad and using blip_to_grad and
  !!   get_matrix_elements.
  !!
  !!   The onsite terms are then replaced by an analytic integration
  !!   which doesn't use the integration grid. This is done by get_onsite_T
  !!
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   4/3/97 CMG
  !!  MODIFICATION HISTORY
  !!   15/08/01 TM
  !!    Added call for get_matrix_elements_new
  !!   29/12/00 TM
  !!    Added call to blip_to_grad_new
  !!   14/05/01 DRB
  !!    Added ROBODoc headers and changed get_matrix_elements calls
  !!   17/05/2001 dave
  !!    Shortened call to blip_to_grad
  !!   17/05/2001 dave
  !!    Shortened variable list in subroutine call
  !!   07/06/2001 dave
  !!    Included in H_matrix_module
  !!   11/06/2001 dave
  !!    Changed to use GenComms for gsum
  !!   09:24, 2003/03/11 dave
  !!    Removed dependence on data_K and calculation of energy, tidied usage
  !!   11:51, 17/03/2003 drb
  !!    Bug fix: added back in data_work
  !!   08:05, 2003/09/22 dave
  !!    Added basis set flag selection; tidied (removed ELEMENTS dependencies)
  !!   16:00, 30/10/2003 drb
  !!    Added PAO assembly call
  !!   2010/03/18 14:27 dave
  !!    Added flag for on-site analytics
  !!   2012/03/13 L.Tong
  !!   - Added spin implementation
  !!   2016/07/13 18:30 nakata
  !!    Renamed H_on_supportfns -> H_on_atomfns
  !!   2016/07/15 18:30 nakata
  !!    Renamed sf_H_sf_rem -> atomf_H_atomf_rem
  !!   2016/07/29 18:30 nakata
  !!    Renamed supports_on_atom -> blips_on_atom
  !!   2016/09/29 18:30 nakata
  !!    Introduced matKEatomf instead of matKE
  !!   2016/12/19 18:30 nakata
  !!    Removed dHrange, matdKE, matdH, paof and sf,
  !!    which are no longer needed.
  !!  SOURCE
  !!
  subroutine get_T_matrix

    use numbers
    use primary_module,              only: bundle
    use matrix_data,                 only: mat, Hrange, aHa_range
    use set_bucket_module,           only: rem_bucket, atomf_H_atomf_rem
    use calc_matrix_elements_module, only: get_matrix_elements_new
    use blip_grid_transform_module,  only: blip_to_grad_new
    use GenComms,                    only: gsum, cq_abort, myid
    use GenBlas,                     only: axpy
    use global_module,               only: flag_basis_set,       &
                                           PAOs, blips,          &
                                           flag_onsite_blip_ana, &
                                           nspin
    use build_PAO_matrices,          only: assemble_2
    use mult_module,                 only: allocate_temp_matrix, &
                                           free_temp_matrix,     &
                                           matrix_sum,           &
                                           matKEatomf
    use functions_on_grid,           only: H_on_atomfns
    use support_spec_format,         only: blips_on_atom
    use species_module,              only: nsf_species

    implicit none

    ! local variables
    integer :: direction, i, np, nn, matwork, this_nsf, spec, spin


    if(flag_basis_set == blips) then
       matwork = allocate_temp_matrix(Hrange,0)
       do direction = 1, 3
          call blip_to_grad_new(myid, direction, H_on_atomfns(1))
          call get_matrix_elements_new(myid, rem_bucket(atomf_H_atomf_rem), &
                                       matwork, H_on_atomfns(1), &
                                       H_on_atomfns(1))
          call matrix_sum(one, matKEatomf, one, matwork)
       end do
       ! replace the onsite blocks with the analytic values...
       if (flag_onsite_blip_ana) then
          i = 1
          do np = 1, bundle%groups_on_node
             if(bundle%nm_nodgroup(np) > 0) then
                do nn = 1, bundle%nm_nodgroup(np)
                   spec = bundle%species(i)
                   this_nsf = nsf_species(spec)
                   call get_onsite_T(blips_on_atom(i), matKEatomf, np, &
                                     nn, i, this_nsf, spec)
                   i = i + 1
                end do
             endif
          end do
       end if
       call free_temp_matrix(matwork)
    else if (flag_basis_set == PAOs) then
       ! Use assemble to generate matrix elements
       call assemble_2(aHa_range, matKEatomf, 2)
    else
       call cq_abort('get_T_matrix: basis set incorrectly specified ', &
                     flag_basis_set)
    end if
    return
  end subroutine get_T_matrix
  !!***

! -----------------------------------------------------------
! Subroutine get_onsite_T
! -----------------------------------------------------------

!!****f* H_matrix_module/get_onsite_T *
!!
!!  NAME
!!   get_onsite_T
!!  USAGE
!!
!!  PURPOSE
!!   This routine evalutes the kinetic energy matrix elements for the block
!!   diagonal. This can be done analytically relatively quickly, because the
!!   support functions are represented on the _same_ blip grid, as they are
!!   associated with a single atom.
!!  INPUTS
!!
!!
!!  USES
!!
!!  AUTHOR
!!   C.M.Goringe
!!  CREATION DATE
!!   02/10/96
!!  MODIFICATION HISTORY
!!   16/05/2001 dave
!!    Converted to F90, added ROBODoc, stripped argument list and
!!    removed NSF=4 dependence
!!   21/06/2001 dave
!!    Included in H_matrix_module
!!   2006/07/18 08:08 dave
!!    Changed to use new blip coefficient storage
!!   2006/09/13 07:53 dave
!!    Added species of atom
!!   2011/11/15 08:02 dave
!!    Changes to blip data
!!   2016/09/29 18:30 nakata
!!    Rename matKE to mat_KE (local name)
!!  SOURCE
!!
  subroutine get_onsite_T(blip_co, mat_KE, np, nn, ip, this_nsf, spec)

    use datatypes
    use numbers
    use GenBlas,             only: axpy, copy, scal, gemm
    use blip,                only: blip_info
    use mult_module,         only: store_matrix_value, scale_matrix_value
    use support_spec_format, only: support_function
    use GenComms,            only: cq_abort
    use global_module,       only: area_ops
    use memory_module,       only: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    ! Shared Variables
    integer :: this_nsf, mat_KE, np, nn, ip, spec

    type(support_function) :: blip_co
    ! Local Variables
    integer, parameter :: MAX_D = 3
    integer :: i1,i2
    real(double) :: temp(this_nsf,this_nsf)

    real(double) ::  FAC(0:MAX_D), D2FAC(0:MAX_D)

    real(double), allocatable, dimension(:) :: work1, work2, work3, work4, work5, work6

    integer :: dx, dy, dz, offset, l, at, nsf1, stat

    allocate(work1(blip_info(spec)%FullArraySize*this_nsf), &
         work2(blip_info(spec)%FullArraySize*this_nsf), &
         work3(blip_info(spec)%FullArraySize*this_nsf), &
         work4(blip_info(spec)%FullArraySize*this_nsf), &
         work5(blip_info(spec)%FullArraySize*this_nsf), &
         work6(blip_info(spec)%FullArraySize*this_nsf), &
         STAT=stat)
    if(stat/=0) call cq_abort("Error allocating arrays for onsite KE elements: ",blip_info(spec)%FullArraySize,this_nsf)
    call reg_alloc_mem(area_ops, 6*this_nsf*blip_info(spec)%FullArraySize,type_dbl)
    ! first, we copy the blip functions for this atom onto a cubic grid;
    ! we make this grid 'too big' in order to have a fast routine below.

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
             if (l.ne.0) then
                at = (((dz+offset)*blip_info(spec)%OneArraySize + &
                     (dy+offset))*blip_info(spec)%OneArraySize + &
                     (dx+offset)) * this_nsf
                do nsf1 = 1,this_nsf
                   work1(nsf1+at) = blip_co%supp_func(nsf1)%coefficients(l)
                enddo
             end if
          end do
       end do
    end do

    ! now, for each direction in turn, we need to apply a number of
    ! 'spreading operations'. Do z first... put blip(z) in 2, and
    ! del2blip(z) in 3

    call copy(blip_info(spec)%FullArraySize*this_nsf,work1,1,work2,1)
    call scal(blip_info(spec)%FullArraySize*this_nsf,FAC(0),work2,1)
    do dz = 1, MAX_D
       offset = dz * blip_info(spec)%OneArraySize * blip_info(spec)%OneArraySize * this_nsf
       call axpy((blip_info(spec)%FullArraySize*this_nsf-offset), FAC(dz), &
            work1(1:), 1, work2(1+offset:), 1 )
       call axpy((blip_info(spec)%FullArraySize*this_nsf-offset), FAC(dz), &
            work1(1+offset:), 1, work2(1:), 1 )
    end do

    call copy(blip_info(spec)%FullArraySize*this_nsf,work1,1,work3,1)
    call scal(blip_info(spec)%FullArraySize*this_nsf,D2FAC(0),work3,1)
    do dz = 1, MAX_D
       offset = dz * blip_info(spec)%OneArraySize * blip_info(spec)%OneArraySize * this_nsf
       call axpy((blip_info(spec)%FullArraySize*this_nsf-offset), D2FAC(dz), &
            work1(1:), 1, work3(1+offset:), 1 )
       call axpy((blip_info(spec)%FullArraySize*this_nsf-offset), D2FAC(dz), &
            work1(1+offset:), 1, work3(1:), 1 )
    end do


    ! now do y : put blip(y).blip(z) in 4,
    ! blip(y).del2blip(z) + del2blip(y).blip(z)  in 5

    call copy(blip_info(spec)%FullArraySize*this_nsf,work2,1,work4,1)
    call scal(blip_info(spec)%FullArraySize*this_nsf,FAC(0),work4,1)
    do dy = 1, MAX_D
       offset = dy * blip_info(spec)%OneArraySize * this_nsf
       call axpy((blip_info(spec)%FullArraySize*this_nsf-offset), FAC(dy), &
            work2(1:), 1, work4(1+offset:), 1 )
       call axpy((blip_info(spec)%FullArraySize*this_nsf-offset), FAC(dy), &
            work2(1+offset:), 1, work4(1:), 1 )
    end do

    call copy(blip_info(spec)%FullArraySize*this_nsf,work2,1,work5,1)
    call scal(blip_info(spec)%FullArraySize*this_nsf,D2FAC(0),work5,1)
    do dy = 1, MAX_D
       offset = dy * blip_info(spec)%OneArraySize * this_nsf
       call axpy((blip_info(spec)%FullArraySize*this_nsf-offset), D2FAC(dy), &
            work2(1:), 1, work5(1+offset:), 1 )
       call axpy((blip_info(spec)%FullArraySize*this_nsf-offset), D2FAC(dy), &
            work2(1+offset:), 1, work5(1:), 1 )
    end do

    call axpy(blip_info(spec)%FullArraySize*this_nsf,FAC(0),work3,1,work5,1)
    do dy = 1, MAX_D
       offset = dy * blip_info(spec)%OneArraySize * this_nsf
       call axpy((blip_info(spec)%FullArraySize*this_nsf-offset), FAC(dy), &
            work3(1:), 1, work5(1+offset:), 1 )
       call axpy((blip_info(spec)%FullArraySize*this_nsf-offset), FAC(dy), &
            work3(1+offset:), 1, work5(1:), 1 )
    end do

    ! and x - put it all into 6

    call copy(blip_info(spec)%FullArraySize*this_nsf,work5,1,work6,1)
    call scal(blip_info(spec)%FullArraySize*this_nsf,FAC(0),work6,1)
    do dx = 1, MAX_D
       offset = dx * this_nsf
       call axpy((blip_info(spec)%FullArraySize*this_nsf-offset), FAC(dx), &
            work5(1:), 1, work6(1+offset:), 1 )
       call axpy((blip_info(spec)%FullArraySize*this_nsf-offset), FAC(dx), &
            work5(1+offset:), 1, work6(1:), 1 )
    end do

    call axpy(blip_info(spec)%FullArraySize*this_nsf,D2FAC(0),work4,1,work6,1)
    do dx = 1, MAX_D
       offset = dx * this_nsf
       call axpy((blip_info(spec)%FullArraySize*this_nsf-offset), D2FAC(dx), &
            work4(1:), 1, work6(1+offset:), 1 )
       call axpy((blip_info(spec)%FullArraySize*this_nsf-offset), D2FAC(dx), &
            work4(1+offset:), 1, work6(1:), 1 )
    end do

    ! and now get the matrix elements by multiplication...
    temp = zero
    call gemm('n','t',this_nsf,this_nsf,blip_info(spec)%OneArraySize*blip_info(spec)%OneArraySize*blip_info(spec)%OneArraySize, &
         one,work1,this_nsf,work6,this_nsf,zero,temp,this_nsf )
    call start_timer(tmr_std_matrices)
    do i1 = 1,this_nsf
       do i2=1,this_nsf
          call scale_matrix_value(mat_KE,np,nn,ip,0,i2,i1,zero,1)
          call store_matrix_value(mat_KE,np,nn,ip,0,i2,i1,blip_info(spec)%SupportGridSpacing*temp(i2,i1),1)
       end do
    end do
    call stop_timer(tmr_std_matrices)
    deallocate(work1,work2,work3, work4,work5,work6, STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating arrays for onsite KE blip elements: ",blip_info(spec)%FullArraySize,this_nsf)
    call reg_dealloc_mem(area_ops, 6*this_nsf*blip_info(spec)%FullArraySize,type_dbl)
    return
  end subroutine get_onsite_T
!!***

  ! -----------------------------------------------------------
  ! Subroutine get_output_energies
  ! -----------------------------------------------------------

  !!****f* H_matrix_module/get_output_energies *
  !!
  !!  NAME
  !!   get_output_energies
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Calculate the Hartree and XC energies for a given density
  !!   Used to find the correct DFT energy which requires the Ha
  !!   and XC energies for the output densities
  !!   Also finds the local PS/neutral atom energy, which integrates
  !!   the potential with the density (as opposed to taking the trace
  !!   of K with the appropriate matrix).
  !!
  !!   It's not immediately obvious where to put this subroutine, but
  !!   it's in H_matrix_module because it closely parallels the energy
  !!   calculations in get_h_on_atomfns
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   2021/07/22
  !!  MODIFICATION HISTORY
  !!   2021/07/26 10:54 dave
  !!    Added local PS for nonNA runs, tidied output
  !!  SOURCE
  !!
  subroutine get_output_energies(rho, size)

    use datatypes
    use numbers
    use global_module,          only: flag_pcc_global, nspin, spin_factor, &
                                      flag_neutral_atom, area_ops, iprint_ops
    use hartree_module,         only: get_hartree_energy
    use XC,                     only: get_xc_energy
    use energy,                 only: hartree_energy_total_rho,  &
                                      xc_energy,       &
                                      local_ps_energy, &
                                      hartree_energy_drho, hartree_energy_drho_input
    use density_module,         only: density_pcc, density_atom
    use maxima_module,          only: maxngrid
    use memory_module,          only: reg_alloc_mem,              &
                                      reg_dealloc_mem, type_dbl
    use GenComms,               only: cq_abort, gsum, inode, ionode
    use pseudopotential_common, only: pseudopotential, flag_neutral_atom_projector
    use dimens,                 only: grid_point_volume, n_my_grid_points
    use GenBlas,                only: dot

    implicit none

    ! Passed variables
    integer :: size
    real(double), dimension(:,:) :: rho

    ! Local variables
    real(double), allocatable, dimension(:,:) :: density_wk
    integer :: spin, stat

    ! Workspace to store various forms of density
    allocate(density_wk(size,nspin), STAT=stat)
    if (stat /= 0) call cq_abort("Error allocating density_wk: ", stat)
    call reg_alloc_mem(area_ops, size * nspin, type_dbl)
    ! Construct total density in density_wk
    density_wk = zero
    do spin = 1, nspin
       density_wk(:,1) = density_wk(:,1) + spin_factor*rho(:,spin)
    end do
    ! Calculate Hartree energy of density (or drho for NA)
    if( flag_neutral_atom ) then
       if(.NOT.flag_neutral_atom_projector) then
          local_ps_energy = &
               dot(n_my_grid_points, pseudopotential, 1, density_wk(:,1), 1) * &
               grid_point_volume
          call gsum(local_ps_energy)
          if(inode==ionode .and. iprint_ops>3) &
               write(io_lun,fmt='(2x,"Output density neutral atom energy (delta rho): ",f19.12)') local_ps_energy
       end if
       hartree_energy_drho_input = hartree_energy_drho
       density_wk(:,1) = density_wk(:,1) - density_atom(:)
       call get_hartree_energy(density_wk, maxngrid, hartree_energy_drho)
       if(inode==ionode .and. iprint_ops>3) &
            write(io_lun,fmt='(2x,"Output density Hartree energy (delta rho): ",f19.12)') hartree_energy_drho
    else
       local_ps_energy = &
            dot(n_my_grid_points, pseudopotential, 1, density_wk(:,1), 1) * &
            grid_point_volume
       call gsum(local_ps_energy)
       if(inode==ionode .and. iprint_ops>3) &
            write(io_lun,fmt='(2x,"Output density local PS energy (delta rho): ",f19.12)') local_ps_energy
       call get_hartree_energy(density_wk, maxngrid, hartree_energy_total_rho)
       if(inode==ionode .and. iprint_ops>3) &
            write(io_lun,fmt='(2x,"Output density Hartree energy (delta rho): ",f19.12)') hartree_energy_drho
    end if
    !  
    ! Construct XC energy of density
    if (flag_pcc_global) then
       density_wk = zero
       do spin = 1, nspin
          density_wk(:,spin) = rho(:,spin) + half * density_pcc(:)
       end do
       call get_xc_energy(density_wk, xc_energy, size)
    else
       call get_xc_energy(rho, xc_energy, size)
    end if
    if(inode==ionode .and. iprint_ops>3) &
         write(io_lun,fmt='(2x,"Output density XC energy: ",f19.12)') xc_energy
    deallocate(density_wk)
    call reg_dealloc_mem(area_ops, size*nspin, type_dbl)
  end subroutine get_output_energies
!!***

  ! -----------------------------------------------------------
  ! Subroutine matrix_scale_diag
  ! -----------------------------------------------------------

  !!****f* H_matrix_module/matrix_scale_diag *
  !!
  !!  NAME
  !!   matrix_scale_diag
  !!  USAGE
  !!
  !!  PURPOSE
  !!   Scales the projector-localised orbital scalar
  !!   product matrix
  !!  INPUTS
  !!
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   D.R.Bowler
  !!  CREATION DATE
  !!   25/6/98
  !!  MODIFICATION HISTORY
  !!   21/06/2001 dave
  !!    Included into H_matrix_module and reduced use list
  !!   24/06/2002 dave
  !!    Added explicit loop over NSF and an if clause for second (zeroing) loop
  !!    Also added comments to end of if and do statements to make more readable
  !!  SOURCE
  !!
  subroutine matrix_scale_diag(matAP, species, n_projectors, l_core, &
                               recip_scale, range)

    use datatypes
    use numbers,        only: zero
    use group_module,   only: parts
    use primary_module, only: bundle
    use cover_module,   only: BCS_parts
    use matrix_data,    only: mat
    use mult_module,    only: scale_matrix_value

    implicit none

    integer :: matAP, species(:), n_projectors(:),l_core( :, : )
    integer :: range

    real(double) :: recip_scale( :, : )

    ! Local variables
    integer :: i,k,n, species_k,l2,np,nb,isu,ind_cover,ind_qart, nsf1, ip

    ! Loop over the i elements and apply the appropriate scaling factor
    ip = 0
    call start_timer(tmr_std_matrices)
    do np = 1,bundle%groups_on_node
       if(bundle%nm_nodgroup(np)>0) then
          do i = 1,bundle%nm_nodgroup(np)
             ip = ip + 1
             if(mat(np,range)%n_nab(i)>0) then
                do nb = 1,mat(np,range)%n_nab(i)
                   isu = mat(np,range)%i_acc(i)+nb-1
                   ind_cover=mat(np,range)%i_part(isu)
                   ind_qart=BCS_parts%lab_cell(ind_cover)
                   k = parts%icell_beg(ind_qart)+mat(np,range)%i_seq(isu)-1
                   species_k = species(k)
                   do n = 1, n_projectors(species_k)
                      l2 = l_core( n, species_k )
                      do nsf1 = 1,mat(np,range)%ndimi(i)
                         call scale_matrix_value(matAP, np, i, ip, nb, &
                                                 nsf1, n,              &
                                                 recip_scale(l2,       &
                                                 species_k))
                      end do ! End do nsf1
                   enddo ! End do n=n_projectors
                enddo ! End do nb=mat%n_nab(i)
             endif ! End if(mat%n_nab(i)>0)
          enddo ! End do i=bundle%nm_nodgroup(np)
       endif ! End if(bundle%nm_nodgroup>0)
    enddo ! End do np=bundle%groups_on_node
    call stop_timer(tmr_std_matrices)
    return
  end subroutine matrix_scale_diag
  !!***

end module H_matrix_module
