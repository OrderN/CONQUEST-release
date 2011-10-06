! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
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
!!   atoms, blip_grid_transform_module, blip, calc_matrix_elements_module, common,
!!   cover_module, datatypes, dimens, generic_blas, generic_comms, grid_index, group_module,
!!   matrix_data, maxima_module, multiply_module, mult_module, numbers, primary_module,
!!   pseudopotential_data, set_blipgrid_module, set_bucket_module, species_module, workspace
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
!!    Removed all explicit references to data_ variables and rewrote in terms of new 
!!    matrix routines
!!   2006/03/06 04:46 dave
!!    Removing passed arrays support, workspace_support
!!   2008/02/01 17:47 dave
!!    Changes for output to file not stdout
!!   2008/05/22 ast
!!    Added timers
!!  TODO
!!   08:28, 2003/09/22 dave
!!    Understand and document onsite_T
!!  SOURCE
!!
module H_matrix_module

  use global_module, ONLY: io_lun
  use timer_stdclocks_module, ONLY: start_timer,stop_timer,tmr_std_hmatrix,tmr_std_allocation,tmr_std_matrices

  implicit none

  ! RCS tag for object file identification
  character(len=80), save, private :: RCSid = "$Id$"
  logical :: locps_output
  integer :: locps_choice

!!***

contains

! ------------------------------------------------------------------------------
! Subroutine get_H_matrix
! ------------------------------------------------------------------------------

!!****f* H_matrix_module/get_H_matrix *
!!
!!  NAME 
!!   get_H_matrix
!!  USAGE
!! 
!!  PURPOSE
!!   gets the entire Hamiltonian matrix
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   C.M.Goringe
!!  CREATION DATE
!!   11/10/95 
!!  MODIFICATION HISTORY
!!   10/08/2000 TM
!!    New get_matrix_elements
!!   14/05/2001 DRB
!!    Added ROBODoc header and removed unnecessary arguments
!!   17/05/2001 dave
!!    Reduced call to get_T_matrix
!!    Reduced call to get_h_on_support
!!    Reduced call to get_HNL_matrix
!!   18/05/2001 dave
!!    Reduced overall argument list
!!   07/06/2001 dave
!!    Included in H_matrix_module
!!   11/06/2001 dave
!!    Changed to use GenComms for gsum
!!   13/05/2002 dave
!!    Added more comments and shifted format statements around
!!   15/11/2002 tm
!!    Changed pseudopotential_data use to pseudopotential_common
!!   09:11, 2003/03/11 dave
!!    Removed find_chdens argument and writing out of energies and
!!    added flag to force or not the rebuilding of data_NL and data_KE
!!  (the non-local and kinetic energy contributions to the H matrix)
!!   2008/05/22 ast
!!    Added timers
!!   2011/10/06 13:51 dave
!!    Added cDFT call to store H without constraint
!!  SOURCE
!!
  subroutine get_H_matrix(rebuild_KE_NL, fixed_potential, electrons, rho, size)

    use datatypes
    use numbers
    use matrix_data, ONLY: dHrange, Hrange, Srange
    use mult_module, ONLY: matNL, matKE, matH, matdH, allocate_temp_matrix, free_temp_matrix, &
         matrix_scale, matrix_sum, matS
    use pseudopotential_common, ONLY: non_local, pseudopotential
    use set_bucket_module, ONLY: rem_bucket, sf_H_sf_rem, pao_H_sf_rem
    use calc_matrix_elements_module, ONLY: get_matrix_elements_new
    use GenComms, ONLY: gsum, end_comms, my_barrier, inode, ionode
    use global_module, ONLY: iprint_ops, flag_vary_basis, flag_basis_set, PAOs, &
         sf, paof, IPRINT_TIME_THRES1,iprint_SC, flag_perform_cDFT
    use PAO_grid_transform_module, ONLY: single_PAO_to_grid
    use functions_on_grid, ONLY: supportfns, H_on_supportfns, &
         allocate_temp_fn_on_grid, free_temp_fn_on_grid, gridfunctions, fn_on_grid
    use io_module, ONLY: dump_matrix, dump_blips, dump_charge, write_matrix, dump_charge
    use potential_module, ONLY: potential
    use dimens, ONLY: n_my_grid_points

    use group_module, ONLY: parts
    use primary_module, ONLY: bundle
    use cover_module, ONLY: BCS_parts
    !use io_module, ONLY: dump_matrix, write_matrix
    use timer_module            ! This is used to declare a local timer
    use cdft_data, ONLY: matHzero,cDFT_NumberAtomGroups,cDFT_Vc,matWc

    implicit none

    ! Passed variables
    logical :: rebuild_KE_NL, fixed_potential

    real(double) :: electrons
    integer :: size
    real(double), dimension(size) :: rho
    
    ! local variables
    real(double) :: kinetic_energy, nl_energy
    integer :: pao_support
    integer :: length, stat, paolength, matwork
    type(cq_timer) :: tmr_l_hmatrix

    call start_timer(tmr_std_hmatrix)                ! Total
    call start_timer(tmr_l_hmatrix,WITH_LEVEL)   ! Just this call

    if(inode==ionode.AND.iprint_ops>3) write(io_lun,fmt='(10x,"Entering get_H_matrix")')
    if(flag_vary_basis.AND.flag_basis_set==PAOs) matwork = allocate_temp_matrix(dHrange,0,paof,sf)
    call matrix_scale(zero,matH)
    if(flag_vary_basis.AND.flag_basis_set==PAOs) call matrix_scale(zero,matdH)
    gridfunctions(H_on_supportfns)%griddata = zero

    if(rebuild_KE_NL) then
       if(inode==ionode.AND.iprint_ops>3) write(io_lun,fmt='(2x,"Rebuilding KE")')
       call matrix_scale(zero,matKE)
       call matrix_scale(zero,matNL)
       ! get the T matrix and the kinetic energy...
       call get_T_matrix(matKE)
       ! now, we do the non-local part (if we are doing it)
       if ( non_local ) then
          if(inode==ionode.AND.iprint_ops>3) write(io_lun,fmt='(2x,"Rebuilding NL")')
          call get_HNL_matrix(matNL)
       end if
       if(iprint_ops>4) call dump_matrix("NNL",matNL,inode)
       if(iprint_ops>4) call dump_matrix("NKE",matKE,inode)
    end if
    ! from here on, workspace support becomes h_on_support...
    ! in fact, what we are getting here is (H_local - T) acting on support
    call get_h_on_support( iprint_ops, fixed_potential, electrons, rho, size)
    if(inode==ionode.AND.iprint_ops>2) write(io_lun,*) 'Doing integration'
    ! Do the integration - support holds <phi| and workspace_support holds H|phi>
    call get_matrix_elements_new(inode-1,rem_bucket(sf_H_sf_rem),matH,supportfns,H_on_supportfns)
    if(inode==ionode.AND.iprint_ops>2) write(io_lun,*) 'Done integration'
    if(iprint_ops>2) call dump_matrix("Nl",matH,inode)
    ! After doing this, we need \chi_{ilm} projected onto the grid, and then a new call to get_matrix_elements_new
    if(flag_vary_basis.AND.flag_basis_set==PAOs) then
       ! Project PAOs onto grid
       if(inode==ionode.AND.iprint_ops>2) write(io_lun,*) 'Doing single_pao_on_support'
       pao_support = allocate_temp_fn_on_grid(paof)
       call single_PAO_to_grid(pao_support)
       ! Do integration
       if(inode==ionode.AND.iprint_ops>2) write(io_lun,*) 'Doing integration'
       call get_matrix_elements_new(inode-1,rem_bucket(pao_H_sf_rem),matwork,pao_support,H_on_supportfns)
       if(inode==ionode.AND.iprint_ops>2) write(io_lun,*) 'Done integration'
       call my_barrier()
       if(inode==ionode.AND.iprint_ops>2) write(io_lun,*) 'Doing axpy'
       call matrix_sum(one,matdH,one,matwork)
       call free_temp_fn_on_grid(pao_support)
       call free_temp_matrix(matwork)
       if(inode==ionode.AND.iprint_ops>2) write(io_lun,*) 'Done axpy'
    end if

    ! add the kinetic energy and non-local matrices to give the complete H matrix
    call matrix_sum(one,matH,half,matKE)
    call matrix_sum(one,matH,one,matNL)
    if(iprint_ops>4) call dump_matrix("NS",matS,inode)
    if(iprint_ops>4) call dump_matrix("NH",matH,inode)
    if(iprint_SC>2) call dump_charge(rho,size,inode)
    ! Store the new H matrix for cDFT
    if(flag_perform_cDFT) then
       call matrix_sum(zero,matHzero,one,matH)
    endif

    call stop_print_timer(tmr_l_hmatrix, "get_H_matrix", IPRINT_TIME_THRES1)
    call stop_timer(tmr_std_hmatrix)

    return
  end subroutine get_H_matrix
!!***

! -----------------------------------------------------------
! Subroutine get_h_on_support
! -----------------------------------------------------------

!!****f* H_matrix_module/get_h_on_support *
!!
!!  NAME 
!!   get_h_on_support
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
!!    Removed find_chdens and total_energy arguments and call to get_electronic_density
!!   2004/10/05 drb
!!    Added hartree_module use
!!   13:58, 2006/07/17
!!    Added functional type selector
!!   2008/04/02  M. Todorovic
!!    Added local potential output
!!   2010/11/10  TM
!!    Moved the parts of the calculation of delta_E_xc from subroutines get_**_xc_potential
!!    to get_h_on_support for the implementation of PCC (partial core correction).
!!   2011/03/31  M.Arita
!!    Added the contribution from P.C.C.
!!  SOURCE
!!
  subroutine get_h_on_support( output_level, fixed_potential, electrons, rho, size)

    use datatypes
    use numbers
    use global_module, ONLY: sf, &
                             flag_functional_type, functional_lda_pz81, &
                             functional_lda_gth96, functional_lda_pw92, &
                             functional_gga_pbe96, functional_gga_pbe96_rev98, &
                             functional_gga_pbe96_r99, area_ops, flag_pcc_global
    use GenBlas, ONLY: copy, axpy, dot, rsum
    use dimens, ONLY: grid_point_volume, n_my_grid_points, n_grid_z
    use block_module, ONLY: n_blocks, n_pts_in_block
    use primary_module, ONLY: domain
    use set_blipgrid_module, ONLY: naba_atm
    use density_module, ONLY : density_pcc
    use GenComms, ONLY: gsum, inode, ionode, cq_abort
    use energy, ONLY: hartree_energy, xc_energy, local_ps_energy, delta_E_hartree
    use hartree_module, ONLY: hartree
    use functions_on_grid, ONLY: gridfunctions, fn_on_grid, H_on_supportfns, supportfns
    use calc_matrix_elements_module, ONLY: norb
    use pseudopotential_common, ONLY: pseudopotential
    use potential_module, ONLY: potential
    use maxima_module, ONLY: maxngrid
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl
    use fft_module, ONLY: fft3, hartree_factor, z_columns_node, i0
    use io_module, ONLY: dump_locps
    use energy, ONLY: delta_E_xc

    implicit none

    ! Shared variables
    logical :: fixed_potential

    integer :: output_level, size, igrid

    real(double) :: electrons
    real(double), dimension(size) :: rho

    ! Local variables
    integer :: n, m, nb, atom, nsf1, point, stat, i, pot_flag
    ! Allocated arrays
    real(double) :: fften
    real(double), allocatable, dimension(:) :: xc_potential, h_potential
    real(double), allocatable, dimension(:) :: xc_epsilon   ! energy_density of XC 
    real(double), allocatable, dimension(:) :: density_wk   ! rho(:)+density_pcc(:)
    logical :: dump_pot(4)

    complex(double_cplx), allocatable, dimension(:) :: chdenr, locpotr

    !call dump_locps(pseudopotential,size,inode)
    !allocate(chdenr(size), locpotr(size), STAT=stat)
    !call fft3( rho, chdenr, size, -1 )
    !call fft3( pseudopotential, locpotr, size, -1 )
    !fften = zero
    !do i = 1, z_columns_node(inode)*n_grid_z
    !   if(i==i0) then
    !      write(io_lun,*) 'G=0: ',grid_point_volume*real(chdenr(i),double)*real(locpotr(i),double), &
    !           aimag(chdenr(i))*aimag(locpotr(i))*grid_point_volume
    !      !locpotr(i0)=cmplx(zero,zero)
    !   else
    !      fften = fften + real(chdenr(i),double)*real(locpotr(i),double)*grid_point_volume - &
    !           aimag(chdenr(i))*aimag(locpotr(i))*grid_point_volume
    !   end if
    !end do
    !write(io_lun,*) 'Energy via FFT: ',fften
    !call fft3( pseudopotential, locpotr, size, 1 )
    !deallocate(chdenr, locpotr, STAT=stat)
    call start_timer(tmr_std_allocation)
    !old allocate(xc_potential(n_my_grid_points), h_potential(maxngrid), STAT=stat)
    allocate(xc_potential(size), xc_epsilon(size), h_potential(size), STAT=stat)
    if(stat/=0) call cq_abort("Error allocating potentials in H: ",n_my_grid_points, stat)
    !oldcall reg_alloc_mem(area_ops, n_my_grid_points+maxngrid,type_dbl)
    call reg_alloc_mem(area_ops, 3*maxngrid,type_dbl)
    call stop_timer(tmr_std_allocation)
    !     first initialise some arrays
    gridfunctions(H_on_supportfns)%griddata = zero
    potential = zero
    h_potential = zero
    xc_potential = zero

    electrons = grid_point_volume * rsum( n_my_grid_points, rho, 1 )
    call gsum(electrons)
    if (inode==ionode.and.output_level>=1) write(io_lun,114) electrons

    !     now calculate the hartree potential on the grid
    call hartree( rho, h_potential, maxngrid, hartree_energy )
    ! Correction term
    delta_E_hartree = -hartree_energy

    ! for P.C.C.
    if (flag_pcc_global) then
      allocate(density_wk(size))
      !density_wk = zero
      density_wk = rho + density_pcc
    endif

    select case(flag_functional_type)
       case (functional_lda_pz81)
          if (flag_pcc_global) then
            call get_xc_potential( density_wk, xc_potential, xc_epsilon, xc_energy, size )
          else
            call get_xc_potential( rho, xc_potential, xc_epsilon, xc_energy, size )
          endif
       case (functional_lda_gth96)
          if (flag_pcc_global) then
            call get_GTH_xc_potential( density_wk, xc_potential, xc_epsilon, xc_energy, size )
          else
            call get_GTH_xc_potential( rho, xc_potential, xc_epsilon, xc_energy, size )
          endif
       case (functional_lda_pw92)
          if (flag_pcc_global) then
            call get_xc_potential_LDA_PW92( density_wk, xc_potential, xc_epsilon, xc_energy, size )
          else
            call get_xc_potential_LDA_PW92( rho, xc_potential, xc_epsilon, xc_energy, size )
          endif
       case (functional_gga_pbe96)
          if (flag_pcc_global) then
            call get_xc_potential_GGA_PBE( density_wk, xc_potential, xc_epsilon, xc_energy, size )
          else
            call get_xc_potential_GGA_PBE( rho, xc_potential, xc_epsilon, xc_energy, size )
          endif
       case (functional_gga_pbe96_rev98)
          if (flag_pcc_global) then
            call get_xc_potential_GGA_PBE( density_wk, xc_potential, xc_epsilon, xc_energy, size, functional_gga_pbe96_rev98 )
          else
            call get_xc_potential_GGA_PBE( rho, xc_potential, xc_epsilon, xc_energy, size, functional_gga_pbe96_rev98 )
          endif
       case (functional_gga_pbe96_r99)
          if (flag_pcc_global) then
            call get_xc_potential_GGA_PBE( density_wk, xc_potential, xc_epsilon, xc_energy, size, functional_gga_pbe96_r99 )
          else
            call get_xc_potential_GGA_PBE( rho, xc_potential, xc_epsilon, xc_energy, size, functional_gga_pbe96_r99 )
          endif
       case default
          !ORI call get_xc_potential( rho, xc_potential, xc_energy, n_my_grid_points )
          if (flag_pcc_global) then
            call get_xc_potential( density_wk, xc_potential, xc_epsilon, xc_energy, size )
          else
            call get_xc_potential( rho, xc_potential, xc_epsilon, xc_energy, size )
          endif
    end select
    ! Calculation of delta_E_xc     ---- 2010.Oct.30 TM 
    !    moved from get_xc_... to here,  by introducing 
    !    xc_epsilon(size) [XC energy density, not potential] 

    ! delta_E_xc 
    delta_E_xc = zero
     if (flag_pcc_global) then
       do igrid = 1, n_my_grid_points
         delta_E_xc = delta_E_xc + density_wk(igrid) * xc_epsilon(igrid) - rho(igrid) * &
                      xc_potential(igrid)
       enddo
     else
       do igrid=1,n_my_grid_points
         delta_E_xc = delta_E_xc + (xc_epsilon(igrid)-xc_potential(igrid))*rho(igrid)
       enddo
     endif
    call gsum(delta_E_xc)
    delta_E_xc = delta_E_xc*grid_point_volume

    ! Make total potential
    if (.not.fixed_potential) then
       call copy( n_my_grid_points, h_potential, 1, potential, 1 )
       call axpy( n_my_grid_points, one, xc_potential, 1, potential, 1 )
       call axpy( n_my_grid_points, one, pseudopotential, 1, potential, 1 )
    end if

    ! Print potential, if necessary
    if (locps_output) then
       dump_pot = (/.false.,.false.,.false.,.false./)
       if (locps_choice<0.OR.locps_choice>15) then           
          if (inode==ionode) write (io_lun,*) 'Bad choice for local potential printout: no output.'    
       else
          if (inode==ionode) write (io_lun,*) 'Writing local potential to file(s).' 
          do i=4,1,-1 
             pot_flag = mod(locps_choice/(2**(i-1)),2)
             if (pot_flag>0) dump_pot(i) = .true. 
             locps_choice = locps_choice - pot_flag*(2**(i-1))
          end do
       end if
       if (dump_pot(1)) call dump_locps("Hartree",h_potential,size,inode)
       if (dump_pot(2)) call dump_locps("XC",xc_potential,size,inode)
       if (dump_pot(3)) call dump_locps("PS",pseudopotential,size,inode)
       if (dump_pot(4)) call dump_locps("Total",potential,size,inode)
    end if

    ! get the pseudopotential energy
    local_ps_energy = dot( n_my_grid_points, pseudopotential, 1, rho, 1 ) * &
         grid_point_volume
    call gsum(local_ps_energy)

    ! now act with the potential on the support functions to get
    ! (H_local - T) on support
    n = 0
    m = 0
    do nb=1,domain%groups_on_node
       if(naba_atm(sf)%no_of_atom(nb)>0) then
          do atom = 1, naba_atm(sf)%no_of_atom(nb)
             do nsf1 = 1, norb(naba_atm(sf),atom,nb)
                do point = 1, n_pts_in_block
                   n = n + 1
                   gridfunctions(H_on_supportfns)%griddata(n) = gridfunctions(supportfns)%griddata(n)*potential(m+point)
                end do
             end do
          end do
       end if ! (naba_atm(sf)%no_of_atom(nb)>0)
       m = m + n_pts_in_block
    end do
114 format(20x,'Electron Count         : ',f30.15)
    deallocate(xc_potential,h_potential, STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating potentials in H: ",n_my_grid_points, stat)
    if (flag_pcc_global) then
      deallocate(density_wk, STAT=stat)
      if (stat .NE. 0) call cq_abort("Error deallocating density_wk: ", n_my_grid_points, stat)
    endif
    !old call reg_dealloc_mem(area_ops, n_my_grid_points+maxngrid,type_dbl)
    call reg_dealloc_mem(area_ops, 3*maxngrid,type_dbl)
    return
  end subroutine get_h_on_support
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
!!    Added lines to allow calculation of derivative of <\phi|\hat{NL}|\phi> with PAOs as basis for \phi
!!  SOURCE
!!
  subroutine get_HNL_matrix(matNL)

    use datatypes
    use numbers
    use matrix_data, ONLY: mat, SPrange, PAOPrange, dHrange
    use mult_module, ONLY: mult, SP_PS_H, PAOP_PS_H, SP_trans, matCS, matdH, matSC, &
         allocate_temp_matrix, free_temp_matrix, matrix_product, matrix_sum, matrix_transpose
    use pseudopotential_data, ONLY: n_projectors, l_core, recip_scale
    use pseudopotential_common, ONLY: pseudo_type, OLDPS, SIESTA, STATE, ABINIT
    use species_module, ONLY: species
    use set_bucket_module, ONLY: rem_bucket, sf_nlpf_rem
    use calc_matrix_elements_module, ONLY: get_matrix_elements_new
    use global_module, ONLY: flag_basis_set, PAOs, blips, flag_vary_basis, paof, nlpf, sf, iprint_ops
    use GenComms, ONLY: cq_abort, myid, inode, ionode
    use GenBlas, ONLY: axpy
    use build_PAO_matrices, ONLY: assemble_2
    use functions_on_grid, ONLY: supportfns, pseudofns
    use io_module, ONLY: dump_matrix

    implicit none

    ! Passed variables
    integer :: matNL

    ! Local variables
    integer :: stat, matdSC, matdCNL, matSCtmp

    if(flag_vary_basis.AND.flag_basis_set==PAOs) then
       matdSC = allocate_temp_matrix(PAOPrange,SP_trans,paof,nlpf)
       matdCNL = allocate_temp_matrix(dHrange,0,paof,sf)
    end if
    matSCtmp = allocate_temp_matrix(SPrange,SP_trans,sf,nlpf)
    ! first, get the overlap of support functions with core pseudowavefunctions
    if(flag_basis_set==blips) then
       if(inode==ionode.AND.iprint_ops>2) write(io_lun,*) 'Calling get_matrix_elements'
       call get_matrix_elements_new(myid,rem_bucket(sf_nlpf_rem),matSC,supportfns,pseudofns)
       !call dump_matrix("NSC",matSC,inode)
    else if(flag_basis_set==PAOs) then
       ! Use assemble to generate matrix elements
       if(inode==ionode.AND.iprint_ops>2) write(io_lun,*) 'Calling assemble'
       if(flag_vary_basis) then
          call assemble_2(SPrange, matSC,3,matdSC)
       else
          call assemble_2(SPrange, matSC,3)       
       end if
       if(inode==ionode.AND.iprint_ops>2) write(io_lun,*) 'Called assemble'
    else
       call cq_abort('get_HNL_matrix: basis set incorrectly specified ',flag_basis_set)
    end if
    if(inode==ionode.AND.iprint_ops>2) write(io_lun,*) 'Made SP'
    call matrix_sum(zero,matSCtmp,one,matSC)
    if(mult(SP_PS_H)%mult_type==2) then ! type 2 means no transpose necessary
       select case(pseudo_type)
       case(OLDPS) 
          call matrix_scale_diag(matSC, species, n_projectors, l_core, recip_scale,SPrange)
       case(SIESTA)
          call matrix_scale_diag_tm(matSC,SPrange)
       case(ABINIT)
          call matrix_scale_diag_tm(matSC,SPrange)
       end select
       call matrix_product(matSCtmp,matSC,matNL,mult(SP_PS_H))
       if(inode==ionode.AND.iprint_ops>2) write(io_lun,*) 'Calling mult_wrap'
       if(flag_vary_basis.AND.flag_basis_set==PAOs) then
          call matrix_product(matdSC,matSC,matdCNL,mult(PAOP_PS_H))
       end if
    else ! Transpose SP->PS, then mult
       if(inode==ionode.AND.iprint_ops>2) write(io_lun,*) 'Type 1 ',matSC,matCS
       call matrix_transpose(matSC,matCS)
       if(inode==ionode.AND.iprint_ops>2) write(io_lun,*) 'Done transpose' 
       select case(pseudo_type)
       case(OLDPS) 
          call matrix_scale_diag(matSC, species, n_projectors, l_core, recip_scale,SPrange)
       case(SIESTA)
          if(inode==ionode.AND.iprint_ops>2) write(io_lun,*) 'Doing scale' 
          call matrix_scale_diag_tm(matSC,SPrange)
          if(inode==ionode.AND.iprint_ops>2) write(io_lun,*) 'Calling scale'
          if(flag_vary_basis.AND.flag_basis_set==PAOs) call matrix_scale_diag_tm(matdSC,PAOPrange)
       case(ABINIT)
          if(inode==ionode.AND.iprint_ops>2) write(io_lun,*) 'Doing scale' 
          call matrix_scale_diag_tm(matSC,SPrange)
          if(inode==ionode.AND.iprint_ops>2) write(io_lun,*) 'Calling scale'
          if(flag_vary_basis.AND.flag_basis_set==PAOs) call matrix_scale_diag_tm(matdSC,PAOPrange)
       end select
       call matrix_product(matSC,matCS,matNL,mult(SP_PS_H))
       if(inode==ionode.AND.iprint_ops>2) write(io_lun,*) 'Calling mult_wrap'
       if(flag_vary_basis.AND.flag_basis_set==PAOs) then 
          call matrix_product(matdSC,matCS,matdCNL,mult(PAOP_PS_H))
       end if
    endif
    if(flag_vary_basis.AND.flag_basis_set==PAOs) then
       if(inode==ionode.AND.iprint_ops>2) write(io_lun,*) 'PAOs'
       call matrix_sum(one,matdH,one,matdCNL)
    end if
    call free_temp_matrix(matSCtmp)
    if(flag_vary_basis.AND.flag_basis_set==PAOs) then
       call free_temp_matrix(matdCNL)
       call free_temp_matrix(matdSC)
    end if
    return

  contains

    ! Added nonef to allow scaling of the matrix <PAO|projector> required for PAO basis of support functions
    ! DRB 2004/07/23 08:28
    subroutine matrix_scale_diag_tm(matSC,range)

      ! Module usage
      use datatypes, ONLY: double
      use numbers, ONLY: very_small, zero
      use group_module, ONLY: parts
      use primary_module, ONLY: bundle
      use cover_module, ONLY: BCS_parts
      use matrix_data, ONLY: mat
      use species_module, ONLY: species
      use pseudo_tm_info, ONLY: pseudo
      use mult_module, ONLY: scale_matrix_value, store_matrix_value
      use GenComms, ONLY: inode, ionode

      implicit none

      ! Passed variables
      integer :: range, matSC

      ! Local variables
      integer :: np, i, nb, isu, ind_cover, ind_qart, ip
      integer :: k, species_k, nl, nsf1
      integer :: mmax, mm, n_proj

      ! Loop over the i elements and apply the appropriate scaling factor
      ip = 0
      call start_timer(tmr_std_matrices)
      do np = 1,bundle%groups_on_node
         if(bundle%nm_nodgroup(np)>0) then
            do i = 1,bundle%nm_nodgroup(np)
               ip = ip+1
               if(mat(np,range)%n_nab(i)>0) then
                  do nb = 1,mat(np,range)%n_nab(i)
                     isu = mat(np,range)%i_acc(i)+nb-1
                     ind_cover=mat(np,range)%i_part(isu)
                     ind_qart=BCS_parts%lab_cell(ind_cover)
                     k = parts%icell_beg(ind_qart)+mat(np,range)%i_seq(isu)-1
                     species_k = species(k)
                     n_proj=0
                    if(pseudo(species_k)%n_pjnl > 0) then
                     do nl = 1, pseudo(species_k)%n_pjnl
                        mmax=pseudo(species_k)%pjnl_l(nl)*2+1
                        do mm=1,mmax
                           n_proj=n_proj+1
                           do nsf1=1,mat(np,range)%ndimi(i)
                              call scale_matrix_value(matSC,np,i,ip,nb,nsf1,n_proj,pseudo(species_k)%pjnl_ekb(nl))
                              if(abs(pseudo(species_k)%pjnl_ekb(nl)) < very_small.AND.inode==ionode) &
                                   write(io_lun,*) 'ekb = 0!!   for nl, species_k, ekb = ', &
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
!!  SOURCE
!!
  subroutine get_T_matrix(matKE)

    use numbers
    use primary_module, ONLY : bundle
    use matrix_data, ONLY : mat, Hrange, dHrange
    use set_bucket_module, ONLY: rem_bucket, sf_H_sf_rem
    use calc_matrix_elements_module, ONLY: get_matrix_elements_new
    use blip_grid_transform_module, ONLY: blip_to_grad_new
    use GenComms, ONLY: gsum, cq_abort, myid
    use GenBlas, ONLY: axpy
    use global_module, ONLY: flag_basis_set, PAOs, blips, flag_vary_basis, paof, sf, flag_onsite_blip_ana
    use build_PAO_matrices, ONLY: assemble_2
    use mult_module, ONLY: allocate_temp_matrix, free_temp_matrix, matdH, matrix_sum
    use functions_on_grid, ONLY: H_on_supportfns
    use support_spec_format, ONLY: supports_on_atom
    use species_module, ONLY: nsf_species

    implicit none

    ! Passed variables
    integer :: matKE

    ! local variables
    integer :: direction, i, np, nn, matwork, matdKE, this_nsf, spec

    matwork = allocate_temp_matrix(Hrange,0)
    if(flag_basis_set==PAOs.AND.flag_vary_basis) matdKE = allocate_temp_matrix(dHrange,0,paof,sf)
    if(flag_basis_set==blips) then
       do direction = 1, 3
          call blip_to_grad_new(myid, direction, H_on_supportfns)
          call get_matrix_elements_new(myid,rem_bucket(sf_H_sf_rem),matwork,H_on_supportfns,H_on_supportfns)
          call matrix_sum(one,matKE,one,matwork)
       end do
       ! replace the onsite blocks with the analytic values...
       if(flag_onsite_blip_ana) then
          i = 1
          do np = 1,bundle%groups_on_node
             if(bundle%nm_nodgroup(np)>0) then
                do nn=1,bundle%nm_nodgroup(np)
                   spec = bundle%species(i)
                   this_nsf = nsf_species(spec)
                   call get_onsite_T(supports_on_atom(i), matKE, np,nn, i,this_nsf,spec)
                   i = i+1
                end do
             endif
          end do
       end if
    else if(flag_basis_set==PAOs) then
       ! Use assemble to generate matrix elements
       if(flag_vary_basis) then
          call assemble_2(Hrange, matKE,2,matdKE)
          call matrix_sum(one,matdH,half,matdKE)
       else
          call assemble_2(Hrange, matKE,2)
       end if
    else
       call cq_abort('get_T_matrix: basis set incorrectly specified ',flag_basis_set)
    end if
    if(flag_basis_set==PAOs.AND.flag_vary_basis) call free_temp_matrix(matdKE)
    call free_temp_matrix(matwork)
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
!!  SOURCE
!!
  subroutine get_onsite_T(blip_co,matKE, np,nn, ip,this_nsf,spec)

    use datatypes
    use numbers
    use GenBlas, ONLY: axpy, copy, scal, gemm
    use blip, ONLY: blip_info, BlipArraySize, OneArraySize, FullArraySize, SupportGridSpacing
    use mult_module, ONLY: store_matrix_value, scale_matrix_value
    use support_spec_format, ONLY: support_function
    use GenComms, ONLY: cq_abort
    use global_module, ONLY: area_ops
    use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl

    implicit none

    ! Shared Variables
    integer :: this_nsf, matKE,np, nn, ip, spec

    type(support_function) :: blip_co
    ! Local Variables
    integer, parameter :: MAX_D = 3
    integer :: i1,i2
    real(double) :: temp(this_nsf,this_nsf)

    real(double) ::  FAC(0:MAX_D), D2FAC(0:MAX_D)

    real(double), allocatable, dimension(:) :: work1, work2, work3, work4, work5, work6

    integer :: dx, dy, dz, offset, l, at, nsf1, stat

    allocate(work1(FullArraySize(spec)*this_nsf),work2(FullArraySize(spec)*this_nsf),work3(FullArraySize(spec)*this_nsf), &
         work4(FullArraySize(spec)*this_nsf),work5(FullArraySize(spec)*this_nsf),work6(FullArraySize(spec)*this_nsf), &
         STAT=stat)
    if(stat/=0) call cq_abort("Error allocating arrays for onsite KE elements: ",FullArraySize(spec),this_nsf)
    call reg_alloc_mem(area_ops, 6*this_nsf*FullArraySize(spec),type_dbl)
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
    offset = BlipArraySize(spec)+1

    do dx = -BlipArraySize(spec), BlipArraySize(spec)
       do dy = -BlipArraySize(spec), BlipArraySize(spec)
          do dz = -BlipArraySize(spec), BlipArraySize(spec)
             l = blip_info(spec)%blip_number(dx,dy,dz)
             if (l.ne.0) then
                at = (((dz+offset)*OneArraySize(spec) + (dy+offset))*OneArraySize(spec) + &
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

    call copy(FullArraySize(spec)*this_nsf,work1,1,work2,1)
    call scal(FullArraySize(spec)*this_nsf,FAC(0),work2,1)
    do dz = 1, MAX_D
       offset = dz * OneArraySize(spec) * OneArraySize(spec) * this_nsf
       call axpy((FullArraySize(spec)*this_nsf-offset), FAC(dz), &
            work1(1:), 1, work2(1+offset:), 1 )
       call axpy((FullArraySize(spec)*this_nsf-offset), FAC(dz), &
            work1(1+offset:), 1, work2(1:), 1 )
    end do

    call copy(FullArraySize(spec)*this_nsf,work1,1,work3,1)
    call scal(FullArraySize(spec)*this_nsf,D2FAC(0),work3,1)
    do dz = 1, MAX_D
       offset = dz * OneArraySize(spec) * OneArraySize(spec) * this_nsf
       call axpy((FullArraySize(spec)*this_nsf-offset), D2FAC(dz), &
            work1(1:), 1, work3(1+offset:), 1 )
       call axpy((FullArraySize(spec)*this_nsf-offset), D2FAC(dz), &
            work1(1+offset:), 1, work3(1:), 1 )
    end do


    ! now do y : put blip(y).blip(z) in 4,
    ! blip(y).del2blip(z) + del2blip(y).blip(z)  in 5

    call copy(FullArraySize(spec)*this_nsf,work2,1,work4,1)
    call scal(FullArraySize(spec)*this_nsf,FAC(0),work4,1)
    do dy = 1, MAX_D
       offset = dy * OneArraySize(spec) * this_nsf
       call axpy((FullArraySize(spec)*this_nsf-offset), FAC(dy), &
            work2(1:), 1, work4(1+offset:), 1 )
       call axpy((FullArraySize(spec)*this_nsf-offset), FAC(dy), &
            work2(1+offset:), 1, work4(1:), 1 )
    end do

    call copy(FullArraySize(spec)*this_nsf,work2,1,work5,1)
    call scal(FullArraySize(spec)*this_nsf,D2FAC(0),work5,1)
    do dy = 1, MAX_D
       offset = dy * OneArraySize(spec) * this_nsf
       call axpy((FullArraySize(spec)*this_nsf-offset), D2FAC(dy), &
            work2(1:), 1, work5(1+offset:), 1 )
       call axpy((FullArraySize(spec)*this_nsf-offset), D2FAC(dy), &
            work2(1+offset:), 1, work5(1:), 1 )
    end do

    call axpy(FullArraySize(spec)*this_nsf,FAC(0),work3,1,work5,1)
    do dy = 1, MAX_D
       offset = dy * OneArraySize(spec) * this_nsf
       call axpy((FullArraySize(spec)*this_nsf-offset), FAC(dy), &
            work3(1:), 1, work5(1+offset:), 1 )
       call axpy((FullArraySize(spec)*this_nsf-offset), FAC(dy), &
            work3(1+offset:), 1, work5(1:), 1 )
    end do

    ! and x - put it all into 6

    call copy(FullArraySize(spec)*this_nsf,work5,1,work6,1)
    call scal(FullArraySize(spec)*this_nsf,FAC(0),work6,1)
    do dx = 1, MAX_D
       offset = dx * this_nsf
       call axpy((FullArraySize(spec)*this_nsf-offset), FAC(dx), &
            work5(1:), 1, work6(1+offset:), 1 )
       call axpy((FullArraySize(spec)*this_nsf-offset), FAC(dx), &
            work5(1+offset:), 1, work6(1:), 1 )
    end do

    call axpy(FullArraySize(spec)*this_nsf,D2FAC(0),work4,1,work6,1)
    do dx = 1, MAX_D
       offset = dx * this_nsf
       call axpy((FullArraySize(spec)*this_nsf-offset), D2FAC(dx), &
            work4(1:), 1, work6(1+offset:), 1 )
       call axpy((FullArraySize(spec)*this_nsf-offset), D2FAC(dx), &
            work4(1+offset:), 1, work6(1:), 1 )
    end do

    ! and now get the matrix elements by multiplication...
    temp = zero
    call gemm('n','t',this_nsf,this_nsf,OneArraySize(spec)*OneArraySize(spec)*OneArraySize(spec), &
         one,work1,this_nsf,work6,this_nsf,zero,temp,this_nsf )
    call start_timer(tmr_std_matrices)
    do i1 = 1,this_nsf
       do i2=1,this_nsf
          call scale_matrix_value(matKE,np,nn,ip,0,i2,i1,zero,1)
          call store_matrix_value(matKE,np,nn,ip,0,i2,i1,SupportGridSpacing(spec)*temp(i2,i1),1)
       end do
    end do
    call stop_timer(tmr_std_matrices)
    deallocate(work1,work2,work3, work4,work5,work6, STAT=stat)
    if(stat/=0) call cq_abort("Error deallocating arrays for onsite KE blip elements: ",FullArraySize(spec),this_nsf)
    call reg_dealloc_mem(area_ops, 6*this_nsf*FullArraySize(spec),type_dbl)
    return
  end subroutine get_onsite_T
!!***

! -----------------------------------------------------------
! Subroutine get_xc_potential
! -----------------------------------------------------------

!!****f* H_matrix_module/get_xc_potential *
!!
!!  NAME
!!   get_xc_potential
!!  USAGE
!!
!!  PURPOSE
!!   Calculates the exchange-correlation potential
!!   on the grid within LDA using the Ceperley-Alder
!!   interpolation formula. It also calculates the
!!   total exchange-correlation energy.
!!
!!   Note that this is the Perdew-Zunger parameterisation of the Ceperley-Alder
!!   results for a homogeneous electron gas, as described in Phys. Rev. B 23, 5048 (1981),
!!   with Ceperley-Alder in Phys. Rev. Lett. 45, 566 (1980)
!!  INPUTS
!!
!!
!!  USES
!!
!!  AUTHOR
!!   E.H.Hernandez
!!  CREATION DATE
!!   02/03/95
!!  MODIFICATION HISTORY
!!   01/11/2005 Antonio
!!    Bibliography:
!!       Exchange energy: It can be found in:
!!                        Phys. Rev. B 45, 13244 (1992)
!!                      **See Eq. 26
!!       Correlation energy: The correlation functional is PZ-81
!!                           Phys. Rev. B 23, 5048 (1981)
!!                         **See Appendix C, and specially Table XII,
!!                             Eq. C1 for the xc hole (rs),
!!                             Eqs. C3 & C4, for rs > 1 and
!!                             Eqs. C5 & C6, for rs < 1
!!   17/05/2001 dave
!!    Converted to F90, added ROBODoc header, shortened call
!!   08/06/2001 dave
!!    Changed to use gsum from GenComms and added RCS Id and Log tags
!!   21/06/2001 dave
!!    Included in H_matrix_module
!!   08:00, 2003/03/12 dave
!!    Added calculation of correction term to band energy
!!   11:49, 30/09/2003 drb
!!    Added comments to make separate testing of X and C easier
!!   2008/03/03 18:32 dave
!!    Removed dsqrt
!!  SOURCE
!!
  subroutine get_xc_potential(density,xc_potential,xc_epsilon, xc_energy,size )

    use datatypes
    use numbers
    use dimens, ONLY: grid_point_volume, one_over_grid_point_volume, &
         n_my_grid_points
    use GenComms, ONLY: gsum

    implicit none

    ! Passed variables
    integer,intent(in) :: size
    real(double),intent(in) :: density(size)
    real(double),intent(out) :: xc_potential(size), xc_epsilon(size)
    real(double),intent(out) :: xc_energy

    !     Local variables
    integer n
    real(double) :: denominator, e_correlation, &
         e_exchange, ln_rs, numerator, &
         rcp_rs, rho, rs, rs_ln_rs, &
         sq_rs, v_correlation, v_exchange
    real(double), parameter :: alpha = -0.45817_double
    real(double), parameter :: beta_1 = 1.0529_double
    real(double), parameter :: beta_2 = 0.3334_double
    real(double), parameter :: gamma = - 0.1423_double
    real(double), parameter :: p = 0.0311_double
    real(double), parameter :: q = - 0.048_double
    real(double), parameter :: r = 0.0020_double
    real(double), parameter :: s = - 0.0116_double
    real(double), parameter :: third = 1.0_double/3.0_double


    xc_energy = zero
    do n=1,n_my_grid_points ! loop over grid pts and store potl on each
       rho = density(n)
       if (rho>very_small) then ! Find radius of hole
          rcp_rs = ( four_thirds * pi * rho )**(third)
       else
          rcp_rs = zero
       end if
       e_exchange = alpha * rcp_rs
       v_exchange = four_thirds * e_exchange
       if (rcp_rs>zero) then
          rs = one/rcp_rs
       else
          rs = zero
       end if
       sq_rs = sqrt(rs)
       if (rs>=one) then
          denominator = one / (one + beta_1 * sq_rs + beta_2 * rs)
          numerator = one + seven_sixths * beta_1 * sq_rs +  &
               four_thirds * beta_2 * rs
          e_correlation = gamma*denominator
          v_correlation = gamma * numerator * denominator * denominator
       else if ((rs<one).and.(rs>very_small)) then
          ln_rs = log(rs)
          rs_ln_rs = rs * ln_rs
          e_correlation = p * ln_rs + q  + r * rs_ln_rs + s * rs
          v_correlation = e_correlation -  &
               third * (p + s * rs + r * (rs_ln_rs + rs))
       else
          e_correlation = zero
          v_correlation = zero
       end if
       ! Both X and C
       xc_energy = xc_energy+(e_exchange+e_correlation)*density(n)
       xc_potential(n) = v_exchange + v_correlation
       xc_epsilon(n)=e_exchange+e_correlation
       ! These two for testing
       ! Just C
       !xc_energy = xc_energy+e_correlation*density(n)
       !xc_potential(n) = v_correlation
       !xc_epsilon(n) = e_correlation
       ! Just X
       !xc_energy = xc_energy+e_exchange*density(n)
       !xc_potential(n) = v_exchange
       !xc_epsilon(n) = e_exchange
    end do ! do n_my_grid_points
    call gsum(xc_energy)
    ! and 'integrate' the energy over the volume of the grid point
    xc_energy = xc_energy * grid_point_volume
    return
  end subroutine get_xc_potential
!!***

! -----------------------------------------------------------
! Subroutine get_GTH_xc_potential
! -----------------------------------------------------------

!!****f* H_matrix_module/get_GTH_xc_potential *
!!
!!  NAME 
!!   get_xc_potential
!!  USAGE
!! 
!!  PURPOSE
!!   Calculates the exchange-correlation potential
!!   on the grid within LDA using the Ceperley-Alder
!!   interpolation formula. It also calculates the
!!   total exchange-correlation energy.
!!
!!   Note that this is the Goedecker/Teter/Hutter formula which 
!!   involves only ratios of polynomials, and is rather easy to 
!!   differentiate.  See PRB 54, 1703 (1996)
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   14:45, 25/03/2003
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine get_GTH_xc_potential(density,xc_potential,xc_epsilon,xc_energy,size )

    use datatypes
    use numbers
    use dimens, ONLY: grid_point_volume, one_over_grid_point_volume, &
         n_my_grid_points
    use GenComms, ONLY: gsum

    implicit none

    ! Passed variables
    integer,intent(in) :: size
    real(double),intent(in) :: density(size)
    real(double),intent(out) :: xc_potential(size),xc_epsilon(size)
    real(double),intent(out) :: xc_energy

    !     Local variables
    integer n
    real(double) :: denominator, e_correlation, &
         e_exchange, ln_rs, numerator, &
         rcp_rs, rho, rs, rs_ln_rs, &
         sq_rs, v_correlation, v_exchange, drs_dRho, t1, t2, dt1, dt2
    real(double), parameter :: a0=0.4581652932831429_double
    real(double), parameter :: a1=2.217058676663745_double
    real(double), parameter :: a2=0.7405551735357053_double
    real(double), parameter :: a3=0.01968227878617998_double
    real(double), parameter :: b1=1.000000000000000_double
    real(double), parameter :: b2=4.504130959426697_double
    real(double), parameter :: b3=1.110667363742916_double
    real(double), parameter :: b4=0.02359291751427506_double
    real(double), parameter :: third = 1.0_double/3.0_double

    xc_energy = zero
    do n=1,n_my_grid_points ! loop over grid pts and store potl on each
       rho = density(n) 
       if (rho>very_small) then ! Find radius of hole
          rcp_rs = ( four*third * pi * rho )**(third)
          rs = one/rcp_rs
          if(rs<0.01_double) write(io_lun,*) 'rs out of range ',n
       else
          rcp_rs = zero
          rs = zero
          write(io_lun,*) 'rho out of range ',n
       end if
       if(rs>zero) then
          drs_dRho = -rs / (3.0 * rho)
          t1 = a0 + rs*(a1 + rs * (a2 + rs * a3))
          t2 = rs * (b1 + rs * (b2 + rs * (b3 + rs * b4)))
          dt1 = a1 + rs * (2.0 * a2 + rs * 3.0 * a3)
          dt2 = b1 + rs * (2.0 * b2 + rs * (3.0 * b3 + rs * 4.0 * b4))
          xc_energy = xc_energy - (t1/t2)*rho
          xc_potential(n) = -(t1/t2) + rho*drs_dRho * (-dt1 / t2 + t1 * dt2 / (t2 * t2))
          xc_epsilon(n) = -t1/t2   ! 2010.Oct.30 TM  
       else
          xc_potential(n) = zero
       end if
    end do ! do n_my_grid_points
    call gsum(xc_energy)
    ! and 'integrate' the energy over the volume of the grid point
    xc_energy = xc_energy * grid_point_volume
    return
  end subroutine get_GTH_xc_potential
!!***


! -----------------------------------------------------------
! Subroutine get_xc_potential_LDA_PW92
! -----------------------------------------------------------

!!****f* H_matrix_module/get_xc_potential_LDA_PW92 *
!!
!!  NAME
!!   get_xc_potential_LDA_PW92
!!  USAGE
!!
!!  PURPOSE
!!   Calculates the exchange-correlation potential
!!   on the grid within LDA using the Ceperley-Alder
!!   interpolation formula. It also calculates the
!!   total exchange-correlation energy.
!!
!!   Note that this is the Perdew-Wang parameterisation of the Ceperley-Alder
!!   results for a homogeneous electron gas, as described in Phys. Rev. B 45, 13244 (1992),
!!   with Ceperley-Alder in Phys. Rev. Lett. 45, 566 (1980)
!!  INPUTS
!!
!!
!!  USES
!!
!!  AUTHOR
!!   A.S. Torralba
!!  CREATION DATE
!!   01/11/05
!!  MODIFICATION HISTORY
!!   2008/03/03 18:32 dave
!!    Removed dsqrt
!!  SOURCE
!!
  subroutine get_xc_potential_LDA_PW92(density,xc_potential,xc_epsilon, xc_energy_total,size)

    use datatypes
    use numbers
    use dimens, ONLY: grid_point_volume, one_over_grid_point_volume, &
         n_my_grid_points
    use GenComms, ONLY: gsum, inode, ionode

    implicit none

    ! Passed variables
    integer,intent(in) :: size
    real(double),intent(in) :: density(size)
    real(double),intent(out) :: xc_epsilon(size), xc_potential(size)
    real(double),intent(out) :: xc_energy_total

    !     Local variables
    integer :: n
    real(double) :: prefactor, postfactor, denominator, &
                    e_correlation, e_exchange, &
                    rcp_rs, rho, rs, sq_rs, &
                    v_correlation, v_exchange, &
                    delta_prefactor, delta_postfactor


    !     From Table I, Phys. Rev. B 45, 13244 (1992), for reference
    real(double), parameter :: alpha  = 1.0421234_double
    real(double), parameter :: alpha1 = 0.21370_double
    real(double), parameter :: beta1  = 7.5957_double
    real(double), parameter :: beta2  = 3.5876_double
    real(double), parameter :: beta3  = 1.6382_double
    real(double), parameter :: beta4  = 0.49294_double
    real(double), parameter :: A      = 0.031091_double

    !     Precalculated constants
    real(double), parameter :: k00 = 1.611991954_double     ! (4*pi/3)**(1/3)
    real(double), parameter :: k01 = -0.458165347_double    ! -3/(2*pi*alpha)
    real(double), parameter :: k02 = -0.062182_double       ! -2*A
    real(double), parameter :: k03 = -0.0132882934_double   ! -2*A*alpha1
    real(double), parameter :: k04 = 0.4723158174_double    ! 2*A*beta1
    real(double), parameter :: k05 = 0.2230841432_double    ! 2*A*beta2
    real(double), parameter :: k06 = 0.1018665524_double    ! 2*A*beta3
    real(double), parameter :: k07 = 0.03065199508_double   ! 2*A*beta4
    real(double), parameter :: k08 = -0.008858862267_double ! 2*k03/3
    real(double), parameter :: k09 = 0.0787193029_double    ! k04/6
    real(double), parameter :: k10 = 0.074361381067_double  ! k05/3
    real(double), parameter :: k11 = 0.0509332762_double    ! k06/2
    real(double), parameter :: k12 = 0.0204346633867_double ! 2*k07/3
    real(double), parameter :: third = 0.333333333_double   ! 1/3

    xc_energy_total = zero
    do n=1,n_my_grid_points ! loop over grid pts and store potl on each
       rho = density(n)

       if (rho > very_small) then ! Find radius of hole
          rcp_rs = k00 * ( rho**third )
       else
          rcp_rs = zero
       end if

       ! ENERGY

       ! Exchange
       e_exchange = k01 * rcp_rs

       ! Correlation
       if (rcp_rs > zero) then
          rs = one/rcp_rs
       else
          rs = zero
       end if
       sq_rs = sqrt(rs)

       prefactor = k02 + k03*rs
       denominator = sq_rs * ( k04 + sq_rs * ( k05 + sq_rs * ( k06 + k07 * sq_rs)))
       if (denominator > zero) then
          postfactor = log( 1 + 1/denominator )
       else
          postfactor = 0
       end if

       e_correlation = prefactor * postfactor

       ! Both exchange and correlation

        xc_epsilon(n) = e_exchange + e_correlation
        xc_energy_total = xc_energy_total + xc_epsilon(n)*rho

       ! POTENTIAL

       ! Exchange

       v_exchange = four_thirds * e_exchange

       ! Correlation
       !   (derivative of rho * e_correlation)

       ! NOTE: delta_prefactor is actually the derivative of rho*prefactor
       !       delta_postfactor is rho times the derivative of postfactor
       delta_prefactor  = k02 + k08*rs
       if (sq_rs > zero) then
          delta_postfactor = sq_rs * ( k09 + sq_rs*(k10 + sq_rs*( k11 + k12 * sq_rs ))) &
                           / ( denominator * ( 1 + denominator ) )
       else
          delta_postfactor = 0
       end if

       v_correlation = delta_prefactor * postfactor + prefactor * delta_postfactor

       xc_potential(n) = v_exchange + v_correlation


    end do ! do n_my_grid_points
    call gsum(xc_energy_total)
    ! and 'integrate' the energy over the volume of the grid point
    xc_energy_total = xc_energy_total * grid_point_volume

    return
  end subroutine get_xc_potential_LDA_PW92
!!***


! -----------------------------------------------------------
! Subroutine get_xc_potential_GGA_PBE
! -----------------------------------------------------------

!!****f* H_matrix_module/get_xc_potential_GGA_PBE *
!!
!!  NAME
!!   get_xc_potential_GGA_PBE
!!  USAGE
!!
!!  PURPOSE
!!   Calculates the exchange-correlation potential
!!   on the grid within GGA using the
!!   Perdew-Burke-Ernzerhof. It also calculates the
!!   total exchange-correlation energy.
!!
!!   Note that this is the functional described in
!!   Phys. Rev. Lett. 77, 3865 (1996)
!!
!!   It is also (depending on an optional parameter)
!!     either revPBE, Phys. Rev. Lett. 80, 890 (1998), 
!!     or RPBE, Phys. Rev. B 59, 7413 (1999)
!!  INPUTS
!!
!!
!!  USES
!!
!!  AUTHOR
!!   A.S. Torralba
!!  CREATION DATE
!!   11/11/05
!!  MODIFICATION HISTORY
!!   17/07/06 rgradient(1,:) used for storage of final result, before transform,
!!            instead of direct transform of the sum 
!!            rgradient(1,:) + rgradient(2,:) + rgradient(3,:)
!!          (this was less efficient)
!!   15:55, 27/04/2007 drb 
!!     Changed recip_vector, grad_density to (n,3) for speed
!!   2008/11/13 ast
!!     Added revPBE and RPBE
!!  SOURCE
!!
  subroutine get_xc_potential_GGA_PBE(density,xc_potential,xc_epsilon, xc_energy,size,flavour)

    use datatypes
    use numbers
    use global_module, ONLY: rcellx, rcelly, rcellz, &
                             functional_gga_pbe96_rev98, functional_gga_pbe96_r99
    use dimens, ONLY: grid_point_volume, one_over_grid_point_volume, &
         n_my_grid_points, n_grid_x, n_grid_y, n_grid_z
    use GenComms, ONLY: gsum, cq_abort, inode, ionode
    use fft_module, ONLY: fft3, recip_vector

    implicit none

    ! Passed variables
    integer,intent(in) :: size
    real(double),intent(in) :: density(size)
    real(double),intent(out) :: xc_potential(size), xc_epsilon(size)
    real(double),intent(out) :: xc_energy
    integer,intent(in), optional :: flavour


    !     Local variables
    integer :: n
    integer :: selector

    real(double),allocatable,dimension(:) :: grad_density
    real(double),allocatable,dimension(:,:) :: grad_density_xyz
    real(double)  :: rho, grad_rho, rho1_3, rho1_6, &
                     ks, s, s2, &
                     t, t2, A, At2, &
                     factor0, factor1, factor2, factor3, factor4, &
                     denominator0, numerator1, denominator1, num_den1, numerator2, &
                     dt2_drho, dA_drho, &
                     dnumerator1_drho, ddenominator1_drho, &
                     dnumerator1_dt2, dnumerator1_dA, &
                     ddenominator1_dt2, ddenominator1_dA, &
                     dfactor2_dt2, dfactor2_dnumerator1, dfactor2_ddenominator1, &
                     dfactor2_drho, dfactor3_dfactor2, dfactor3_drho
    real(double) :: xc_energy_lda_total, &
                    e_correlation_lda, &
                    de_correlation_lda, &
                    e_exchange, e_correlation, &
                    de_exchange, de_correlation, &
                    dde_exchange, dde_correlation
    real(double) :: df_dgrad_rho
    real(double) :: kappa, mu_kappa
    real(double) :: rpbe_exp

    complex(double_cplx), allocatable, dimension(:,:) :: rgradient      ! Gradient in reciprocal space

    !     From Phys. Rev. Lett. 77, 3865 (1996)
    real(double), parameter :: mu = 0.21951_double
    real(double), parameter :: beta = 0.066725_double
    real(double), parameter :: gamma = 0.031091_double
    real(double), parameter :: kappa_ori = 0.804_double

    !     From Phys. Rev. Lett. 80, 890 (1998)
    real(double), parameter :: kappa_alt = 1.245_double

    !     Precalculated constants
    real(double), parameter :: mu_kappa_ori = 0.27302_double     ! mu/kappa_ori
    real(double), parameter :: mu_kappa_alt = 0.17631_double     ! mu/kappa_alt
    real(double), parameter :: two_mu = 0.43902_double           ! 2*mu
    real(double), parameter :: beta_gamma = 2.146119_double      ! beta/gamma
    real(double), parameter :: beta_X_gamma = 0.002074546_double ! beta*gamma
    real(double), parameter :: k01 = 0.16162045967_double        ! 1/(2*(3*pi*pi)**(1/3))
    real(double), parameter :: k02 = -0.16212105381_double       ! -3*mu*((4*pi/3)**(1/3))/(2*pi*alpha)
                                                                 ! =mu*k00*k01(LDA_PW92)=mu*k04
    real(double), parameter :: k03 = 1.98468639_double           ! ((4/pi)*(3*pi*pi)**(1/3))**(1/2)
    real(double), parameter :: k04 = -0.738558852965_double      ! -3*((4*pi/3)**(1/3))/(2*pi*alpha) = k00*k01 in LDA_PW92
    real(double), parameter :: k05 = 0.05240415_double           ! -2*k01*k02
    real(double), parameter :: k06 = -0.593801317784_double      ! k04*kappa_ori
    real(double), parameter :: k07 = -0.984745137287_double      ! 4*k04/3
    real(double), parameter :: third = 0.333333333_double        ! 1/3
    real(double), parameter :: seven_thirds = 2.333333333_double ! 7/3

    integer :: stat
    !      Selector options
    integer, parameter :: fx_original    = 1                     ! Used in PBE and revPBE
    integer, parameter :: fx_alternative = 2                     ! Used in RPBE

    ! Choose between PBE or revPBE parameters
    if(PRESENT(flavour)) then
      if(flavour==functional_gga_pbe96_rev98) then
        kappa=kappa_alt
        mu_kappa=mu_kappa_alt
      else
        kappa=kappa_ori
        mu_kappa=mu_kappa_ori
      end if 
    else
      kappa=kappa_ori
      mu_kappa=mu_kappa_ori
    end if

    allocate(grad_density(size), grad_density_xyz(size,3),rgradient(size,3),STAT = stat)
    !initialisation  2010.Oct.30 TM
     grad_density(:)=zero; grad_density_xyz(:,:)=zero; rgradient(:,:) = zero
     xc_epsilon(:) = zero; xc_potential(:) = zero
    if(stat/=0) call cq_abort("Error allocating arrays for PBE functional: ",stat)
    ! Choose functional form
    if(PRESENT(flavour)) then
      if(flavour==functional_gga_pbe96_r99) then
        selector=fx_alternative
      else
        selector=fx_original
      end if
    else
      selector=fx_original
    end if

    ! Build the gradient of the density
    call build_gradient (density, grad_density, grad_density_xyz, size)

    ! Get the LDA part of the functional
    call get_xc_potential_LDA_PW92(density, xc_potential, xc_epsilon, xc_energy_lda_total, size)

    xc_energy = zero
    do n=1,n_my_grid_points ! loop over grid pts and store potl on each
       rho = density(n)
       grad_rho = grad_density(n)

       !!!!!!   XC GGA ENERGY

       ! Exchange

       if (rho > very_small) then
          rho1_3 = rho ** third
          rho1_6 = sqrt (rho1_3)
          s = k01 * grad_rho / (rho ** four_thirds)
          s2 = s * s
          if(selector == fx_alternative) then           ! RPBE
            rpbe_exp = exp(-mu_kappa * s2)
            e_exchange = k06*rho1_3*(1.0_double-rpbe_exp)
          else                                          ! PBE, revPBE
            denominator0 = 1.0 / (1.0 + mu_kappa * s2)
            factor0 = k02 * rho1_3
            factor1 = s2 * denominator0
            ! NOTE: This doesn't look like in Phys. Rev. Lett. 77:18, 3865 (1996)
            !       because the 1 in Fx, has been multiplied by Ex-LDA and is implicit
            !       in xc_energy_lda(n), in the total energy below
            e_exchange = factor0 * factor1
          end if
       else
          e_exchange = zero
       end if

       ! Correlation

       if (rho > very_small) then
          e_correlation_lda = xc_epsilon(n) - k04 * rho1_3

          ! t=grad_rho/(2*rho*ks); ks=sqrt(4*kf/pi); kf=(3*pi*pi*rho)**(1/3); s=grad_rho/(2*rho*kf)
          ks = k03 * rho1_6
          t = grad_rho / ( 2 * ks * rho )
          t2 = t * t

          A = exp(-e_correlation_lda / gamma) - 1.0
          if (A > very_small) then
             A = beta_gamma / A
          else
             A = beta_gamma * BIG
          end if

          At2 = A * t2
          numerator1 = 1.0 + At2
          denominator1 = 1.0 + At2 + At2 * At2
          num_den1 = numerator1 / denominator1

          factor2 = t2 * num_den1
          factor3 = gamma * log(one + beta_gamma * factor2)

          e_correlation = factor3;  !gamma * log( 1.0 + beta_gamma * t2 * num_den1 )
       else
          e_correlation = zero
       end if


       !*ast* TEST-POINT 1
       ! Both exchange and correlation
       xc_energy = xc_energy + (xc_epsilon(n) + e_exchange + e_correlation)*rho
       ! Only LDA part
       !xc_energy = xc_energy + (xc_epsilon(n))*rho
       ! LDA + exchange
       !xc_energy = xc_energy + (xc_epsilon(n) + e_exchange)*rho
       ! Only exchange
       !xc_energy = xc_energy + (e_exchange)*rho
       ! LDA + correlation
       !xc_energy = xc_energy + (xc_epsilon(n) + e_correlation)*rho
       ! Only correlation
       !xc_energy = xc_energy + (e_correlation)*rho

       !!!!!!   POTENTIAL

       !!!   Terms due to df/drho

       ! Exchange

       if (rho > very_small) then
          if(selector == fx_alternative) then           ! RPBE
            de_exchange = k07 * rho1_3 * (kappa - (two_mu * s2 + kappa) * rpbe_exp)
          else                                          ! PBE, revPBE
            de_exchange = four_thirds * factor0 * factor1 * ( 1 - 2* denominator0 )
          end if
       else
          de_exchange = zero
       end if

       ! Correlation

       if (rho > very_small) then
          de_correlation_lda = ( xc_potential(n) &
                               - four_thirds * k04 * rho1_3 &
                               - e_correlation_lda ) / rho

          dt2_drho = -seven_thirds * t2 / rho
          dA_drho  = A * A * exp(-e_correlation_lda / gamma ) * de_correlation_lda / beta
          dnumerator1_dt2   = A
          dnumerator1_dA    = t2
          factor4 = 1.0 + 2 * At2
          ddenominator1_dt2 = A * factor4
          ddenominator1_dA  = t2 * factor4
          dnumerator1_drho   = dnumerator1_dt2 * dt2_drho &
                             + dnumerator1_dA * dA_drho
          ddenominator1_drho = ddenominator1_dt2 * dt2_drho &
                             + ddenominator1_dA * dA_drho
          dfactor2_dt2 = num_den1
          dfactor2_dnumerator1   = t2 / denominator1
          dfactor2_ddenominator1 = -factor2 / denominator1
          dfactor2_drho = dfactor2_dt2 * dt2_drho &
                        + dfactor2_dnumerator1 * dnumerator1_drho &
                        + dfactor2_ddenominator1 * ddenominator1_drho
          dfactor3_dfactor2 = beta / ( 1.0 + beta_gamma * factor2 )
          dfactor3_drho = dfactor3_dfactor2 * dfactor2_drho
          de_correlation = factor3 + rho * dfactor3_drho
       else
          de_correlation = zero
       end if

       !*ast* TEST-POINT 2
       xc_potential(n) = xc_potential(n) + de_exchange + de_correlation
       ! Only LDA part
       !xc_potential(n) = xc_potential(n)
       ! LDA + exchange
       !xc_potential(n) = xc_potential(n) + de_exchange
       ! Only exchange
       !xc_potential(n) = de_exchange
       ! LDA + correlation
       !xc_potential(n) = xc_potential(n) + de_correlation
       ! Only correlation
       !xc_potential(n) = de_correlation

       !!!   Terms due to df/d|grad_rho|

       ! Exchange

       if (rho > very_small) then
          if(selector == fx_alternative) then           ! RPBE
             dde_exchange = -k05 * s * rpbe_exp
          else                                          ! PBE, revPBE
             dde_exchange = -k05 * s * denominator0 * denominator0!factor1 * denominator0
          end if 
       else
          dde_exchange = zero
       end if

       ! Correlation

       if (rho > very_small) then
          numerator2 = beta_X_gamma * t * ( 1.0 + 2.0 * At2 ) &
                     / (( gamma * denominator1 + beta * t2 * numerator1 ) * denominator1)
          dde_correlation = numerator2 / ks;
       else
          dde_correlation = zero
       end if

       ! Normalisation (modulus of gradient)

       !*ast* TEST-POINT 3
       if(abs(grad_density(n)) > very_small) then  !DEBUG
       df_dgrad_rho = (dde_exchange + dde_correlation) / grad_density(n)
       else                 !DEBUG
       df_dgrad_rho = zero   !DEBUG
       endif                !DEBUG
       ! Only LDA part
       !df_dgrad_rho = 0.0
       ! LDA + exchange
       !df_dgrad_rho = dde_exchange / grad_density(n)
       ! Only exchange
       !df_dgrad_rho = dde_exchange / grad_density(n)
       ! LDA + correlation
       !df_dgrad_rho = dde_correlation / grad_density(n)
       ! Only correlation
       !df_dgrad_rho = dde_correlation / grad_density(n)


       ! Gradient times derivative of energy

       grad_density_xyz(n,1) = grad_density_xyz(n,1)*df_dgrad_rho
       grad_density_xyz(n,2) = grad_density_xyz(n,2)*df_dgrad_rho
       grad_density_xyz(n,3) = grad_density_xyz(n,3)*df_dgrad_rho

       !xc_epsilon  !2010.Oct.30 TM
       xc_epsilon(n) = xc_epsilon(n) + e_exchange + e_correlation
       !*ast* TEST-POINT 4
       !TM delta_E_xc = delta_E_xc + (xc_epsilon(n) + e_exchange + e_correlation)*rho
       ! Only LDA part
       !delta_E_xc = delta_E_xc + (xc_epsilon(n))*rho
       ! LDA + exchange
       !delta_E_xc = delta_E_xc + (xc_epsilon(n) + e_exchange)*rho
       ! Only exchange
       !delta_E_xc = delta_E_xc + (e_exchange)*rho
       ! LDA + correlation
       !delta_E_xc = delta_E_xc + (xc_epsilon(n) + e_correlation)*rho
       ! Only correlation
       !delta_E_xc = delta_E_xc + (e_correlation)*rho

    end do ! do n_my_grid_points


    !!!   Final steps of the energy calculation

    ! Add the energies and 'integrate' over the volume of the grid points

    call gsum(xc_energy)
    xc_energy = xc_energy * grid_point_volume

    ! Fourier transform the gradient, component by component

    call fft3(grad_density_xyz(:,1), rgradient(:,1), size, -1 )
    call fft3(grad_density_xyz(:,2), rgradient(:,2), size, -1 )
    call fft3(grad_density_xyz(:,3), rgradient(:,3), size, -1 )

    !!!   Get the second term of the potential by taking derivatives

    ! First, get the scalar product of the (normalised) gradient and the wave vector (times i)

    rgradient(:,1) = -rgradient(:,1)*minus_i*recip_vector(:,1)
    rgradient(:,1) = rgradient(:,1) - rgradient(:,2)*minus_i*recip_vector(:,2)
    rgradient(:,1) = rgradient(:,1) - rgradient(:,3)*minus_i*recip_vector(:,3)

    ! Add terms of the scalar product and then Fourier transform back the resultant vector
    ! NOTE: Store the result in the modulus of the gradient (not needed anymore) to save memory

    call fft3( grad_density, rgradient(:,1), size, 1 )

    ! Finally, get the potential
    ! NOTE that here grad_density is NOT the modulus of the gradient,
    !      but a term of the potential (see Fourier transform above)

    do n=1,n_my_grid_points
        xc_potential(n) = xc_potential(n) - grad_density(n)
    end do

    ! I changed the order of deallocation, because I was told that deallocation 
    ! of the latest allocated array should be done first. 
    ! (though I am not sure whether it is true or not)   2010.Oct.30 TM
     
    if(allocated(rgradient)) then
       deallocate(rgradient,STAT=stat)
       if(stat/=0) call cq_abort("Error deallocating rgradient",stat)
    end if
    if(allocated(grad_density_xyz)) then
       deallocate(grad_density_xyz,STAT=stat)
       if(stat/=0) call cq_abort("Error deallocating grad_density_xyz",stat)
    end if
    if(allocated(grad_density)) then
       deallocate(grad_density,STAT=stat)
       if(stat/=0) call cq_abort("Error deallocating grad_density",stat)
    end if

    return
  end subroutine get_xc_potential_GGA_PBE
!!***


! -----------------------------------------------------------
! Subroutine build_gradient
! -----------------------------------------------------------

!!****f* H_matrix_module/build_gradient *
!!
!!  NAME
!!   build_gradient
!!  USAGE
!!
!!  PURPOSE
!!   Calculates the gradient of the density, and its modulus
!!  INPUTS
!!
!!
!!  USES
!!
!!  AUTHOR
!!   A.S. Torralba
!!  CREATION DATE
!!   11/11/05
!!  MODIFICATION HISTORY
!!   15:55, 27/04/2007 drb 
!!    Changed recip_vector, grad_density to (n,3) for speed

!!
  subroutine build_gradient (density, grad_density, grad_density_xyz, size)

    use datatypes
    use numbers
    use global_module, ONLY: rcellx, rcelly, rcellz
    use dimens, ONLY: n_my_grid_points, n_grid_x, n_grid_y, n_grid_z
    use fft_module, ONLY: fft3, recip_vector
    use GenComms, ONLY: cq_abort

    implicit none

    ! Passed variables
    integer,intent(in) :: size
    real(double), intent(in), dimension(size) :: density
    real(double), intent(out), dimension(size) :: grad_density
    !ORI real(double), intent(out), dimension(3, size) :: grad_density_xyz
    real(double), intent(out), dimension(size,3) :: grad_density_xyz

    ! Local variables
    complex(double_cplx), allocatable :: rdensity(:)      ! Density in reciprocal space
    complex(double_cplx), allocatable :: rdensity_tmp(:)  ! Temporal reciprocal density
    integer :: stat

    !local recip_vec_tm
     
    allocate(rdensity(size), rdensity_tmp(size), STAT=stat)
    if(stat /= 0) call cq_abort('ERROR in build_gradient : stat,size = ',stat,size)
    grad_density_xyz = zero               ! to remove SIGFPE   2010.Oct.25 TM
    rdensity = zero ; rdensity_tmp = zero ! to remove SIGFPE   2010.Oct.25 TM

    ! Fourier transform the density
    call fft3( density, rdensity, size, -1 )


    ! Compute the derivative with respect to x
    rdensity_tmp = -rdensity*minus_i*recip_vector(:,1)!*rcellx/n_grid_x
    call fft3( grad_density_xyz(:,1), rdensity_tmp, size, 1 )

    ! Compute the derivative with respect to y
    rdensity_tmp = -rdensity*minus_i*recip_vector(:,2)!*rcelly/n_grid_y
    call fft3( grad_density_xyz(:,2), rdensity_tmp, size, 1 )

    ! Compute the derivative with respect to z
    rdensity_tmp = -rdensity*minus_i*recip_vector(:,3)!*rcellz/n_grid_z
    call fft3( grad_density_xyz(:,3), rdensity_tmp, size, 1 )

    ! Calculate the modulus of the gradient
    grad_density = sqrt(grad_density_xyz(:,1)**2 + grad_density_xyz(:,2)**2 + grad_density_xyz(:,3)**2)


    return
  end subroutine build_gradient
!!***


! -----------------------------------------------------------
! Subroutine map_density
! -----------------------------------------------------------

!!****f* H_matrix_module/map_density *
!!
!!  NAME
!!   map_density
!!  USAGE
!!
!!  PURPOSE
!!   Map cartesian coordinates to density index
!! (Only for one processor and old arrangement of the density)
!! (x,y,z) starts at (1,1,1)
!!  INPUTS
!!
!!
!!  USES
!!
!!  AUTHOR
!!   A.S. Torralba
!!  CREATION DATE
!!   21/11/05
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine map_density (x, y, z, n)

    use block_module, ONLY: nx_in_block, ny_in_block, nz_in_block
    use group_module, ONLY: blocks

    implicit none

    ! Passed variables
    integer, intent(in) :: x, y, z
    integer, intent(out) :: n

    ! Local variables
    integer :: gx, gy, gz, bx, by, bz, n_block

    if((x > (nx_in_block * blocks%ngcellx)).OR.&
       (y > (ny_in_block * blocks%ngcelly)).OR.&
       (z > (nz_in_block * blocks%ngcellz))) &
    then
        n = -1
        return
    end if

    gx=1+mod(x-1, nx_in_block)
    gy=1+mod(y-1, ny_in_block)
    gz=1+mod(z-1, nz_in_block)

    bx=1+(x-1)/nx_in_block
    by=1+(y-1)/ny_in_block
    bz=1+(z-1)/nz_in_block

    if(mod(by+bz, 2) == 1) then
        bx = blocks%ngcellx - bx + 1
    end if

    if(mod(bz, 2) == 0) then
        by = blocks%ngcelly - by + 1
    end if

    n_block = bx + blocks%ngcellx*((by-1) + blocks%ngcelly*(bz-1))

    n = nx_in_block*ny_in_block*nz_in_block*(n_block-1) + &
        gx + nx_in_block*((gy-1) + ny_in_block*(gz-1))

    return
  end subroutine map_density
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
  subroutine matrix_scale_diag(matSC, species, n_projectors, l_core, recip_scale, range)

    use datatypes
    use numbers, ONLY: zero
    use group_module, ONLY: parts
    use primary_module, ONLY: bundle
    use cover_module, ONLY: BCS_parts
    use matrix_data, ONLY : mat
    use mult_module, ONLY: scale_matrix_value

    implicit none

    integer :: matSC, species(:), n_projectors(:),l_core( :, : )
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
                         call scale_matrix_value(matSC,np,i,ip,nb,nsf1,n,recip_scale( l2, species_k ))
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
