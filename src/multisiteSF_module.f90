! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module multisiteSF
! ------------------------------------------------------------------------------
! Code area 11: basis operations
! ------------------------------------------------------------------------------

!!****h* Conquest/multisiteSF *
!!  NAME
!!   multisiteSF
!!
!!  PURPOSE
!!   Contains subroutines related to multi-site support functions.
!!
!!  SUBROUTINES
!!   --- for both SSSFs and MSSFs ---
!!   initial_SFcoeff, normalise_SFcoeff,
!!   --- for on-site SFs (single-site SFs) ---
!!   initial_SFcoeff_onsite, loop_initial_SFcoeff_onsite
!!   --- for Multi-site SFs (MSSFs) ---
!!   smear_MSSF, print_naba_in_MSSF,
!!   --- for Localized filter diagonalisation (LFD) ---
!!   LocFilterDiag,    LFD_make_Subspace_halo, LFD_pickup_subspace_elements, LFD_make_Subspace_i,
!!   LFD_symm_sub_mat, LFD_make_TVEC,          LFD_filter,                   LFD_put_TVEC_to_SFcoeff,
!!   transpose_2Dmat,  LFD_debug_matrix,
!!   LFD_SCF
!!
!!  USES
!!
!!  AUTHOR
!!   A.Nakata
!!
!!  CREATION DATE
!!   2016/10/04
!!
!!  MODIFICATION HISTORY
!!   2019/12/30 tsuyoshi
!!      introduced n_dump_SFcoeff for dumping SFcoeff & (K or L) matrices during the SCF iteration of SFcoeff
!!
!!  SOURCE
!!
module multisiteSF_module

  use datatypes
  use numbers
  use global_module,          only: io_lun, iprint_basis, area_basis
  use GenComms,               only: myid, my_barrier, cq_abort, inode, ionode
  use timer_module,           only: start_timer,stop_timer
  use timer_stdclocks_module, only: tmr_std_basis,tmr_std_matrices
  use memory_module,          only: reg_alloc_mem, reg_dealloc_mem, &
                                    type_int, type_dbl



  implicit none

  ! Number of MSSFs
  logical :: flag_MSSF_nonminimal              !nonmin_mssf
  real(double) :: MSSF_nonminimal_offset       !nonmin_mssf
  integer, allocatable, dimension(:)  :: MSSF_nonminimal_species ! 1:SZ, 2:SZP, 3:others

  ! Smearing MSSF
  logical :: flag_MSSF_smear                                     ! smear MSSF (default=F)
  integer :: MSSF_Smear_type                                     ! smearing-function type (default=1:Fermi-Dirac)
  real(double) :: MSSF_Smear_center
  real(double) :: MSSF_Smear_shift
  real(double) :: MSSF_Smear_width

  ! LFD
  real(double) :: LFD_kT                                         ! kT in LFD filter function (default=0.1)
  real(double) :: LFD_ChemP                                      ! chemical potential in LFD filter function (default=0.0)
  logical :: flag_LFD_useChemPsub                                ! use ChemP of subspaces (default=T)
  logical :: flag_LFD_ReadTVEC                                   ! read trial vectors from input file (default=F)
  real(double), allocatable, dimension(:,:,:) :: LFD_TVEC_read   ! (species, nsf, npao)

  ! LFD SFcoeff iteration
  logical :: flag_LFD_nonSCF                                     ! Update Hpao with SCF density
  logical :: flag_mix_LFD_SCF                                    ! Perform LFD at each SCF cycle
  real(double) :: LFD_threshE                                    ! energy  threshold for SFcoeff iteration by LFD with SCF charge
  real(double) :: LFD_threshD                                    ! density threshold for SFcoeff iteration by LFD with SCF charge
  real(double) :: LFD_Thresh_EnergyRise                          ! energy  threshold when energy rises
  integer :: LFD_max_iteration                                   ! max. iteration of LFD SFcoeff iteration

  ! MD
  logical :: flag_LFD_MD_UseAtomicDensity                        ! use atomic density when recalculate SFcoeff by LFD (default=T)

  ! Frequency for dumping matrices
  integer :: n_dumpSFcoeff
!!***

contains

!!****f* multisiteSF_module/initial_SFcoeff *
!!
!!  NAME
!!   initial_SFcoeff
!!
!!  PURPOSE
!!   Construct initial single/multi-site SF coefficients.
!!   
!!   If flag_LFD is .true., calculate matSpao and matHpao (matKEpao and matNLpao)
!!   and make initial coefficients by sub:LocFilterDiag.
!!
!!   If not, initial coefficients of on-site PAOs are set by calling sub:initial_SFcoeff_onsite
!!
!!   This subroutine is called in sub:initial_phis in initialisation_module.f90
!!   and in sub:update_H in move_atoms.f90.
!!
!!  AUTHOR
!!   A.Nakata
!!  CREATION DATE
!!   14/10/2016
!!  MODIFICATION HISTORY
!!   10/11/2017 nakata
!!    Removed unused grab_charge
!!   02/12/2019 nakata
!!    Removed dump_matrix(SFcoeff), which will be changed to dump_pos_and_matrices in near future
!!  SOURCE
!!
  subroutine initial_SFcoeff(LFD_build_Spao, LFD_build_Hpao, fixed_potential, output_naba_in_MSSF)

    use matrix_data,               only: SFcoeff_range
    use mult_module,               only: mult, matrix_scale, matrix_transpose, matrix_sum, &
                                         matSFcoeff, matSFcoeff_tran, aLa_aHa_aLHa
    use global_module,             only: nspin, spin_factor, nspin_SF,                     &
                                         flag_Multisite, flag_LFD
    use density_module,            only: get_electronic_density, density
    use dimens,                    only: n_my_grid_points
    use maxima_module,             only: maxngrid
    use functions_on_grid,         only: atomfns, H_on_atomfns
    use PAO_grid_transform_module, only: single_PAO_to_grid
    use S_matrix_module,           only: get_S_matrix
    use H_matrix_module,           only: get_H_matrix
    use store_matrix,              only: dump_pos_and_matrices
!    use io_module,                 only: dump_matrix
    use GenComms,                  only: mtime

    implicit none

    ! Passed variables
    logical :: LFD_build_Spao, LFD_build_Hpao, fixed_potential, output_naba_in_MSSF
    ! Local
    real(double), dimension(nspin) :: electrons
    real(double) :: electrons_tot, t0, t1, t2
    integer :: spin_SF


<<<<<<< HEAD
    if (inode==ionode .and. iprint_basis>3) write(io_lun,*) 'We are in sub:initial_SFcoeff'
=======
    if (inode==ionode .and. iprint_basis>1) write(io_lun,*) 'We are in sub:initial_SFcoeff'
>>>>>>> Tidy LFD and start LibXC reference output

    if (output_naba_in_MSSF) call print_naba_in_MSSF

    call my_barrier()
    t0 = mtime()

    if (flag_LFD .and. .not.flag_Multisite) &
       call cq_abort("The LFD method is available only for multi-site SFs.")

    ! Make SF-PAO coefficients
    if (flag_LFD) then
       ! Prepare Spao and Hpao
       ! 1. Calculate PAO values on grids
       if (LFD_build_Spao) call single_PAO_to_grid(atomfns)
       ! 2. Construct matSpao
       if (LFD_build_Spao) call get_S_matrix(inode, ionode, transform_AtomF_to_SF=.false.)
       ! 3. Construct matHpao
       if (LFD_build_Hpao) call get_H_matrix(.true., fixed_potential, electrons, &
                                             density, maxngrid, transform_AtomF_to_SF=.false.)
       call my_barrier()
       t1 = mtime()
<<<<<<< HEAD
       if (inode == ionode .and. iprint_basis>4) write (io_lun,'(A,f20.8,A)') &
=======
       if (inode == ionode .and. iprint_basis>2) write (io_lun,'(A,f20.8,A)') &
>>>>>>> Tidy LFD and start LibXC reference output
          'LFD: Time for S and H in primitive PAO: ', t1-t0, ' ms'
       t0 = t1
       ! Rayson's Localised Filter Diagonalisation method
       call LocFilterDiag(mult(aLa_aHa_aLHa)%ahalo)    
!    else if (flag_orth_MC) then
!       call orth_SFcoeff
    else
       ! Set coefficients (1, 0.1 or 0) only onsite
       ! (available for both conventional single-site and multi-site SFs)
       call initial_SFcoeff_onsite
    endif

    ! Smear multi-site SFcoeff if needed
    if (flag_MSSF_smear) then
       do spin_SF = 1, nspin_SF
          call smear_MSSF(matSFcoeff(spin_SF),SFcoeff_range, &
                          MSSF_Smear_type, MSSF_Smear_center, MSSF_Smear_shift, MSSF_Smear_width)
       enddo
    endif

    ! Normalise
    call normalise_SFcoeff

    ! Transpose
    do spin_SF = 1,nspin_SF
       call matrix_scale(zero,matSFcoeff_tran(spin_SF))
       call matrix_transpose(matSFcoeff(spin_SF), matSFcoeff_tran(spin_SF))
    enddo

    ! Write out current SF coefficients with some iprint (in future)
     if (iprint_basis>=3) call dump_pos_and_matrices(index=98)

    call my_barrier()
    t2 = mtime()
<<<<<<< HEAD
    if (inode == ionode .and. iprint_basis>2) write (io_lun,'(A,f20.8,A)') &
=======
    if (inode == ionode .and. iprint_basis>1) write (io_lun,'(A,f20.8,A)') &
>>>>>>> Tidy LFD and start LibXC reference output
       'Time for making initial multi-site SF coeffficients: ', t2-t0, ' ms'

    return
  end subroutine initial_SFcoeff
!!***

!!****f* multisiteSF_module/normalise_SFcoeff *
!!
!!  NAME 
!!   normalise_SFcoeff
!! 
!!  PURPOSE
!!   Normalise SF coefficients by 1 / sum_b_(Ca,b*Ca,b).
!!
!!   This subroutine is called in sub:
!!
!!  INPUTS
!! 
!!  USES
!! 
!!  AUTHOR
!!   A.Nakata
!!  CREATION DATE
!!   2016/10/27
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine normalise_SFcoeff
    
    use datatypes
    use numbers
    use global_module,  ONLY: id_glob, species_glob, IPRINT_TIME_THRES2, nspin_SF
    use group_module,   ONLY: parts
    use primary_module, ONLY: bundle
    use cover_module,   ONLY: BCS_parts
    use matrix_data,    ONLY: mat, halo, SFcoeff_range
    use mult_module,    ONLY: matSFcoeff, matrix_pos, mat_p
    use pao_format
    use timer_module
    
    implicit none

    ! Local variables
    integer :: iprim, part, memb, neigh, ist
    integer :: gcspart, neigh_global_part, neigh_global_num, spec_j, j_in_halo, count_pao_j
    integer :: spin_SF, nsf_i, sf1, l2, nacz2, m2, wheremat
    real(double) :: val, summ
    type(cq_timer) :: tmr_l_tmp1   


    call start_timer(tmr_std_basis)
    call start_timer(tmr_l_tmp1,WITH_LEVEL)
    if(iprint_basis>=5.AND.myid==0) write(io_lun,'(6x,A)') ' We are in normalise_SFcoeff'

    call start_timer(tmr_std_matrices)

    do spin_SF = 1, nspin_SF
       iprim = 0
       do part = 1,bundle%groups_on_node ! Loop over primary set partitions
          if(bundle%nm_nodgroup(part)>0) then ! If there are atoms in partition
             do memb = 1,bundle%nm_nodgroup(part) ! Loop over atom i
                iprim  = iprim+1
                nsf_i  = mat(part,SFcoeff_range)%ndimi(memb)   ! number of SFs of i
                ! Loop over SFs on i
                do sf1 = 1, nsf_i                 
                   ! Sum up SFcoeff of this SF
                   summ = zero
                   do neigh = 1, mat(part,SFcoeff_range)%n_nab(memb) ! Loop over neighbours j of i
                      ist = mat(part,SFcoeff_range)%i_acc(memb)+neigh-1
                      gcspart = BCS_parts%icover_ibeg(mat(part,SFcoeff_range)%i_part(ist))+mat(part,SFcoeff_range)%i_seq(ist)-1
                      neigh_global_part = BCS_parts%lab_cell(mat(part,SFcoeff_range)%i_part(ist))
                      neigh_global_num  = id_glob(parts%icell_beg(neigh_global_part)+mat(part,SFcoeff_range)%i_seq(ist)-1)
                      spec_j    = species_glob(neigh_global_num)
                      j_in_halo = halo(SFcoeff_range)%i_halo(gcspart)
                      ! Loop over PAOs on j
                      count_pao_j = 1
                      do l2 = 0, pao(spec_j)%greatest_angmom
                         do nacz2 = 1, pao(spec_j)%angmom(l2)%n_zeta_in_angmom
                            do m2 = -l2,l2
                               wheremat = matrix_pos(matSFcoeff(spin_SF),iprim,j_in_halo,sf1,count_pao_j)
                               val = mat_p(matSFcoeff(spin_SF))%matrix(wheremat)
                               summ = summ + val*val
                               count_pao_j = count_pao_j + 1
                            enddo ! m2
                         enddo ! nacz2
                      enddo !l2
                   end do ! neigh
!
                   ! Normalise SFcoeff
                   summ = sqrt(summ)
                   summ = one / summ
                   do neigh = 1, mat(part,SFcoeff_range)%n_nab(memb) ! Loop over j
                      ist = mat(part,SFcoeff_range)%i_acc(memb)+neigh-1
                      gcspart = BCS_parts%icover_ibeg(mat(part,SFcoeff_range)%i_part(ist))+mat(part,SFcoeff_range)%i_seq(ist)-1
                      neigh_global_part = BCS_parts%lab_cell(mat(part,SFcoeff_range)%i_part(ist))
                      neigh_global_num  = id_glob(parts%icell_beg(neigh_global_part)+mat(part,SFcoeff_range)%i_seq(ist)-1)
                      spec_j    = species_glob(neigh_global_num)
                      j_in_halo = halo(SFcoeff_range)%i_halo(gcspart)
                      ! Loop over PAOs on j
                      count_pao_j = 1
                      do l2 = 0, pao(spec_j)%greatest_angmom
                         do nacz2 = 1, pao(spec_j)%angmom(l2)%n_zeta_in_angmom
                            do m2 = -l2,l2
                               wheremat = matrix_pos(matSFcoeff(spin_SF),iprim,j_in_halo,sf1,count_pao_j)
                               val = mat_p(matSFcoeff(spin_SF))%matrix(wheremat)
                               val = val * summ
                               mat_p(matSFcoeff(spin_SF))%matrix(wheremat) = val
                               count_pao_j = count_pao_j + 1
                            enddo ! m2
                         enddo ! nacz2
                      enddo !l2
                   end do ! neigh
                   ! finish this SF
                end do ! sf1
             end do ! memb
          end if ! nm_nodgroup > 0
       end do ! part
    end do ! spin_SF
!
    call stop_timer(tmr_std_matrices)
    call stop_print_timer(tmr_l_tmp1,"normalise_SFcoeff",IPRINT_TIME_THRES2)
    call stop_timer(tmr_std_basis)
!
  return
  end subroutine normalise_SFcoeff
!!***

!!****f* multisiteSF_module/initial_SFcoeff_onsite *
!!
!!  NAME 
!!   initial_SFcoeff_onsite
!! 
!!  PURPOSE
!!   Construct initial single/multi-site SF coefficients.
!!   Initial coefficients of the on-site PAOs are set to some values
!!   and those of the neighbour atoms are set to zero.
!!
!!   Set initial SF coefficients of conventional single-site SFs,
!!   in which primary atom i = neighbour atom j.
!!
!!   This subroutine is based on sub:get_support_pao_rep in ol_rad_table_subs.f90.
!!
!!   For one_to_one SFs (nsf = npao), 
!!   the initial coefficients are set to 1 or 0.
!!   For contracted SFs (nsf /= npao),
!!   the initial coefficients are set to the values which are 
!!   close to 1 for the first zeta, and close to 0 for the other zeta.
!!
!!   This subroutine is called in sub:initial_phis in initialisation_module.f90.
!!
!!  INPUTS
!! 
!!  USES
!! 
!!  AUTHOR
!!   A.Nakata
!!  CREATION DATE
!!   2016/10/04
!!  MODIFICATION HISTORY
!!   2019/12/02 nakata
!!    Removed dump_matrix(SFcoeff), which will be changed to dump_pos_and_matrices in near future
!!
!!  SOURCE
!!
  subroutine initial_SFcoeff_onsite
    
    use datatypes
    use numbers
    use global_module,  ONLY: id_glob, species_glob, IPRINT_TIME_THRES2, nspin_SF
    use store_matrix,   ONLY: dump_pos_and_matrices
!    use io_module,      ONLY: dump_matrix
    use group_module,   ONLY: parts
    use primary_module, ONLY: bundle
    use cover_module,   ONLY: BCS_parts
    use matrix_data,    ONLY: mat, halo, SFcoeff_range
    use mult_module,    ONLY: matSFcoeff, matSFcoeff_tran, matrix_scale, matrix_transpose
    use timer_module
    
    implicit none

    ! Local variables
    integer :: part, memb, neigh, ist, atom_num, atom_spec, iprim
    integer :: gcspart, neigh_global_part, neigh_global_num, neigh_species
    real(double) :: dx, dy, dz, r2
    type(cq_timer) :: tmr_l_tmp1   
    integer :: spin_SF


    call start_timer(tmr_std_basis)
    call start_timer(tmr_l_tmp1,WITH_LEVEL)
    if(iprint_basis>=5.AND.myid==0) write(io_lun,'(6x,A)') ' We are in initial_SFcoeff_onsite'


    if(iprint_basis>=5.AND.myid==0) write(io_lun,'(6x,i5,A)') myid, ' Zeroing matSFcoeff'
    do spin_SF=1,nspin_SF
       call matrix_scale(zero,matSFcoeff(spin_SF))
       call matrix_scale(zero,matSFcoeff_tran(spin_SF))
    enddo
    if(iprint_basis>=5.AND.myid==0) write(io_lun,'(6x,A)') ' Done Zeroing'

    iprim = 0
    call start_timer(tmr_std_matrices)
    do part = 1,bundle%groups_on_node ! Loop over primary set partitions
       if(iprint_basis>=6.AND.myid==0) write(io_lun,fmt='(6x,"Processor, partition: ",2i7)') myid,part
       if(bundle%nm_nodgroup(part)>0) then ! If there are atoms in partition
          do memb = 1,bundle%nm_nodgroup(part) ! Loop over atoms
             atom_num = bundle%nm_nodbeg(part)+memb-1
             iprim=iprim+1
             ! Atomic species
             atom_spec = bundle%species(atom_num)
             if(iprint_basis>=6.AND.myid==0) write(io_lun,'(6x,"Processor, atom, spec: ",3i7)') myid,memb,atom_spec
             do neigh = 1, mat(part,SFcoeff_range)%n_nab(memb) ! Loop over neighbours of atom
                ist = mat(part,SFcoeff_range)%i_acc(memb)+neigh-1
                ! Build the distances between atoms - needed for phases 
                gcspart = BCS_parts%icover_ibeg(mat(part,SFcoeff_range)%i_part(ist))+mat(part,SFcoeff_range)%i_seq(ist)-1
                ! Displacement vector
                dx = BCS_parts%xcover(gcspart)-bundle%xprim(atom_num)
                dy = BCS_parts%ycover(gcspart)-bundle%yprim(atom_num)
                dz = BCS_parts%zcover(gcspart)-bundle%zprim(atom_num)
                r2 = dx*dx + dy*dy + dz*dz
                r2 = sqrt(r2)
                if(iprint_basis>=6.AND.myid==0) write(80+myid,*) 'dx,y,z,r: ',gcspart,dx,dy,dz,r2
                if (r2.le.RD_ERR) then    ! only on-site
                   ! We need to know the species of neighbour
                   neigh_global_part = BCS_parts%lab_cell(mat(part,SFcoeff_range)%i_part(ist)) 
                   neigh_global_num  = id_glob(parts%icell_beg(neigh_global_part)+mat(part,SFcoeff_range)%i_seq(ist)-1)
                   neigh_species = species_glob(neigh_global_num)
                   if (atom_spec.ne.neigh_species) call cq_abort("Error! Species of atoms i and j are different.")
                   ! Now loop over support functions and PAOs and call routine
                   call loop_initial_SFcoeff_onsite(iprim,halo(SFcoeff_range)%i_halo(gcspart),atom_spec)
                end if ! r2 = 0
             end do ! neigh
          end do ! memb
       end if ! nm_nodgroup > 0
    end do ! part

    ! normalise SFcoeff
    call normalise_SFcoeff

    ! transpose SFcoeff
    do spin_SF = 1,nspin_SF
       call matrix_transpose(matSFcoeff(spin_SF), matSFcoeff_tran(spin_SF))
    enddo

    ! Write out current SF coefficients with some iprint (in future)
    ! if (iprint_basis>=3) call dump_pos_and_matrices
!    if (nspin_SF == 1) then
!       call dump_matrix("SFcoeff",    matSFcoeff(1), inode)
!    else
!       call dump_matrix("SFcoeff_up", matSFcoeff(1), inode)
!       call dump_matrix("SFcoeff_dn", matSFcoeff(2), inode)
!    end if

    call my_barrier

    call stop_timer(tmr_std_matrices)
    call stop_print_timer(tmr_l_tmp1,"initial_SFcoeff_onsite",IPRINT_TIME_THRES2)
    call stop_timer(tmr_std_basis)
!
  return
  end subroutine initial_SFcoeff_onsite
!!***

!!****f* multisiteSF_module/loop_initial_SFcoeff_onsite
!!
!!  NAME
!!   loop_initial_SFcoeff_onsite
!!
!!  PURPOSE
!!   Set initial SF coefficients (SF_atom_i, PAO_atom_j)
!!   for conventional single-site SFs, in which i = j.
!!
!!   This subroutine is used instead of 
!!   sub:get_support_pao_rep in ol_rad_table_subs.f90.
!!
!!   This subroutine is called in sub:initial_SFcoeff_onsite.
!!  
!!  AUTHOR
!!   A.Nakata
!!  CREATION DATE
!!   2016/10/04
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine loop_initial_SFcoeff_onsite(iprim,j_in_halo,atom_spec)

    use global_module,  ONLY: nspin_SF
    use mult_module,    ONLY: matSFcoeff, store_matrix_value_pos, matrix_pos
    use pao_format
    use species_module, ONLY: nsf_species, npao_species

    implicit none

    ! Passed variables
    integer, intent(in) ::  iprim, j_in_halo, atom_spec
    ! Local variables
    integer :: wheremat
    real(double) :: val
    integer :: nsf_i, sf1, l2, ll2, nacz2, m2, spin_SF
    integer :: count, count_pao_j


    nsf_i = nsf_species(atom_spec)
    ! Loop over SFs on i
    if (nsf_species(atom_spec).eq.npao_species(atom_spec)) then
    ! not-contracted SFs
       if (iprint_basis>=6) write(io_lun,*) 'SFs of species', atom_spec, ' are not contracted'
       do sf1 = 1, nsf_i
          count_pao_j = 1
          ! Loop over PAOs on j = i
          do l2 = 0, pao(atom_spec)%greatest_angmom
             do nacz2 = 1, pao(atom_spec)%angmom(l2)%n_zeta_in_angmom
                do m2 = -l2,l2
                   if (sf1==count_pao_j) then
                      val = one
                      do spin_SF = 1, nspin_SF
                         wheremat = matrix_pos(matSFcoeff(spin_SF),iprim,j_in_halo,sf1,count_pao_j)
                         call store_matrix_value_pos(matSFcoeff(spin_SF),wheremat,val)
                      enddo
                   endif
                   count_pao_j = count_pao_j + 1
                enddo ! m2
             enddo ! nacz2
          enddo !l2
       enddo ! sf1
    else
    ! contracted SFs
       if (iprint_basis>=6) write(io_lun,*) 'SFs of species', atom_spec, ' are contracted'
       do sf1 = 1, nsf_i
          count_pao_j = 1 ! which PAO
          ! Loop over PAOs on j = i
          do l2 = 0, pao(atom_spec)%greatest_angmom
             do nacz2 = 1, pao(atom_spec)%angmom(l2)%n_zeta_in_angmom
                ! Count indexes which angular momentum channel we're in
                count = 1
                if (l2>0.AND.pao(atom_spec)%angmom(l2)%n_zeta_in_angmom>0) then
                   do ll2 = 0,l2-1
                      if (pao(atom_spec)%angmom(ll2)%n_zeta_in_angmom>0) count = count + (2*ll2+1)
                   end do
                end if
                do m2 = -l2,l2
                   if (sf1==count) then
                      if (nacz2==1) then
                         val = one
                      else
                         val = 0.1_double
                      endif
                      do spin_SF = 1, nspin_SF
                         wheremat = matrix_pos(matSFcoeff(spin_SF),iprim,j_in_halo,sf1,count_pao_j)
                         call store_matrix_value_pos(matSFcoeff(spin_SF),wheremat,val)
                      enddo
                   endif
                   count = count + 1
                   count_pao_j = count_pao_j + 1
                enddo ! m2
             enddo ! nacz2
          enddo !l2
       enddo ! sf1
    endif
!
  return
  end subroutine loop_initial_SFcoeff_onsite
!!***

!!****f* multisiteSF_module/smear_MSSF *
!!
!!  NAME
!!   smear_MSSF
!!  USAGE
!!
!!  PURPOSE
!!   Smear multi-site SF coefficients depending on distances.
!!
!!   This subroutine is mainly for MSSF matSFcoeff, 
!!   but can be used any kinds of matrices.
!!
!!   (r_center+r_shift) is the center position of the smearing function.
!!   In default, r_center is set to r_MS.
!!   WIDTH determines the width of the smearing function (like kT).
!!
!!   Function to be used for smearing is specified by itype:
!!      itype = 1 ... Fermi-Dirac function (with center+shift and width)
!!              2 ... Error function (with center+shift)
!!
!!   This subroutine is called in sub:initial_SFcoeff.
!!
!!  AUTHOR
!!   A.Nakata
!!  CREATION DATE
!!   2017/01/08
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine smear_MSSF(matA,Arange,itype,r_center,r_shift,WIDTH)

    use datatypes
    use global_module, ONLY: id_glob, species_glob, sf, paof, nlpf
    use species_module, ONLY: species, nsf_species, npao_species
    use pseudo_tm_info, ONLY: pseudo
    use GenComms, ONLY: myid, cq_abort, gcopy

    use matrix_module, ONLY: matrix, matrix_halo
    use group_module, ONLY: parts
    use primary_module, ONLY: bundle
    use cover_module, ONLY: BCS_parts
    use matrix_data, ONLY: mat, halo, rcut, Srange
    use functions, ONLY: erfc_cq
    use mult_module, ONLY: matrix_pos, mat_p

    implicit none

    ! Local variables
    integer :: iprim, np, i, j, ist, atom_num, atom_i, atom_spec
    integer :: neigh_global_part, atom_j, neigh_species, j_in_halo
    integer :: sf1, sf2, nsf1, nsf2, k, l
    integer :: gcspart
    integer :: wheremat
    real(double) :: dx, dy, dz, r2
    real(double) :: val

    real(double) :: width1,EXFRM,FLTR

    ! Passed variables
    integer :: matA, Arange, itype
    real(double) :: r_center, r_shift, WIDTH


    if (iprint_basis>=5.and.inode==ionode) write(io_lun,*) 'We are in sub:smear_MSSF'

    sf1 = mat_p(matA)%sf1_type
    sf2 = mat_p(matA)%sf2_type

!   --- smearing ---

    if (itype.eq.1) width1 = one / width

    if (r_center.eq.zero) r_center = rcut(Srange) * half ! = r_SF = r_pao + r_MS

    iprim = 0; atom_i = 0
       do np = 1,bundle%groups_on_node ! Loop over primary set partitions
          if(bundle%nm_nodgroup(np)>0) then ! If there are atoms in partition
             do i = 1,bundle%nm_nodgroup(np) ! Loop over atoms
                atom_num = bundle%nm_nodbeg(np) + i - 1
                iprim = iprim + 1
                atom_i = bundle%ig_prim(iprim)   ! global number of i
                ! Atomic species
                atom_spec = bundle%species(atom_num)
                if (sf1.eq.sf) then
                   nsf1 = nsf_species(atom_spec)
                else if (sf1.eq.paof) then
                   nsf1 = npao_species(atom_spec)
                else if (sf1.eq.nlpf) then
                   nsf1 = 0
                   if (pseudo(atom_spec)%n_pjnl .gt. 0) then
                      do l = 1, pseudo(atom_spec)%n_pjnl
                         nsf1 = nsf1 + pseudo(atom_spec)%pjnl_l(l)*2 + 1
                      enddo
                   endif
                else
                   call cq_abort("Index should be SF, PAO or NLPF in sub: smear_MSSF")
                endif
                do j = 1, mat(np,Arange)%n_nab(i) ! Loop over neighbours of atom
                   ist = mat(np,Arange)%i_acc(i) + j - 1
                   ! Build the distances between atoms
                   gcspart = BCS_parts%icover_ibeg(mat(np,Arange)%i_part(ist))+mat(np,Arange)%i_seq(ist)-1
                   ! Displacement vector
                   dx = BCS_parts%xcover(gcspart)-bundle%xprim(atom_num)
                   dy = BCS_parts%ycover(gcspart)-bundle%yprim(atom_num)
                   dz = BCS_parts%zcover(gcspart)-bundle%zprim(atom_num)
                   r2 = dx*dx + dy*dy + dz*dz
                   r2 = sqrt(r2)
                   ! check the species of neighbour
                   neigh_global_part = BCS_parts%lab_cell(mat(np,Arange)%i_part(ist))
                   atom_j = id_glob(parts%icell_beg(neigh_global_part)+mat(np,Arange)%i_seq(ist)-1)   ! global number of j
                   j_in_halo = halo(Arange)%i_halo(gcspart)
                   neigh_species = species_glob(atom_j)
                   if (sf2.eq.sf) then
                      nsf2 = nsf_species(neigh_species)
                   else if (sf2.eq.paof) then
                      nsf2 = npao_species(neigh_species)
                   else if (sf2.eq.nlpf) then
                      nsf2 = 0
                      if (pseudo(neigh_species)%n_pjnl .gt. 0) then
                         do l = 1, pseudo(neigh_species)%n_pjnl
                            nsf2 = nsf2 + pseudo(neigh_species)%pjnl_l(l)*2 + 1
                         enddo
                      endif
                   else
                      call cq_abort("Index should be SF, PAO or NLPF in sub: smear_MSSF")
                   endif

                   if (itype.eq.1) then
                      EXFRM = ( r2 - (r_center + r_shift) ) * width1     ! (r - r_mc) / width
                      FLTR  = one / (exp(EXFRM) + one)                   ! f(r) = 1 / ( exp[(r - (r_center+r_shift)) / kT] + 1 )
                   else if (itype.eq.2) then
                      FLTR = half * erfc_cq(r2 - (r_center + r_shift))      ! f(r) = 0.5 * erfc(r - (r_center+r_shift))
                   endif
                   if (iprint_basis>=6) write(io_lun,'(2(A,F10.5))') 'r=',r2,'  FLTR =',FLTR

                   do k = 1, nsf1
                      do l = 1, nsf2
                         wheremat = matrix_pos(matA,iprim,j_in_halo,k,l)
                         val = mat_p(matA)%matrix(wheremat)
                         val = val * FLTR
                         mat_p(matA)%matrix(wheremat) = val
                      enddo ! l
                   enddo ! k
                end do ! j
             end do ! i
          end if ! End if nm_nodgroup > 0
       end do ! np
!
    return
  end subroutine smear_MSSF
!!***

!!****f* multisiteSF_module/print_naba_in_MSSF *
!!
!!  NAME
!!   print_naba_in_MSSF
!!
!!  PURPOSE
!!   Count the number of the neighbor atoms in the multi-site range
!!   in new atomic position and write out.
!!
!!   This subroutine is called in sub:initial_SFcoeff.
!!
!!  AUTHOR
!!   A.Nakata
!!  CREATION DATE
!!   2017/02/20
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine print_naba_in_MSSF

    use primary_module,         only: bundle
    use matrix_data,            only: mat, SFcoeff_range
    use global_module,          only: ni_in_cell
    use GenComms,               only: gsum

    implicit none

    ! Local
    integer :: stat, iprim, atom_i, np, i
    integer, allocatable, dimension(:) :: num_naba_MS   ! storage for number of naba atoms


    allocate(num_naba_MS(ni_in_cell),STAT=stat)
    if(stat/=0) call cq_abort('print_naba_in_MSSF: error allocating num_naba_MS')

    num_naba_MS(:) = 0

    iprim = 0; atom_i = 0  
    do np = 1,bundle%groups_on_node           ! Loop over primary set partitions
       if(bundle%nm_nodgroup(np)>0) then      ! If there are atoms in partition
          do i = 1,bundle%nm_nodgroup(np)     ! Loop over atom_i
             iprim = iprim + 1
             atom_i = bundle%ig_prim(iprim)   ! global number of i
             num_naba_MS(atom_i) = mat(np,SFcoeff_range)%n_nab(i)
          enddo ! i
       endif ! endif of (bundle%nm_nodgroup(np)>0)
    enddo ! np

    ! print out the number of neighbour atoms for each target atom
    call my_barrier()
    call gsum(num_naba_MS,ni_in_cell)
    if (inode==ionode .and. iprint_basis>1) then
       write(io_lun,*) ' --- Number of neighbour atoms in the multi-site range ---'
       write(io_lun,*) '     atom ID      number of neighbour atoms'
       do i = 1, ni_in_cell
          write(io_lun,'(5X,I8,5X,I8)') i, num_naba_MS(i)
       enddo
    endif

    ! deallocate subspace matrices their labels
    deallocate(num_naba_MS,STAT=stat)
    if(stat/=0) call cq_abort('print_naba_in_MSSF: error deallocating num_naba_MS')
!
    return
  end subroutine print_naba_in_MSSF
!!***

!!****f* multisiteSF_module/LocFilterDiag *
!!
!!  NAME
!!   LocFilterDiag
!!
!!  PURPOSE
!!   Construct initial multi-site SF coefficients 
!!   by Rayson's localized filter diagonalisation (LFD) method.
!!        M. J. Rayson and P. R. Briddon, Phys. Rev. B 80, 205104 (2009).
!!        M. J. Rayson, Comput. Phys. Comm. 181, 1051-1056 (2010).
!!
!!   Procedure is as follows:
!!   (1) Subspaces of Hpao and Spao are made for NODEs 
!!       (Hsub and Ssub have dimension(pao_halo,pao_halo)).
!!
!!   (2) For each primary ATOMs,
!!       correponding matrix elements are picked up from Hsub and Ssub and
!!       put into Hsub_i and Ssub_i
!!       (Hsub_i and Ssub_i have dimension (pao_neighbours,pao_neigbours).
!!
!!   (3) Diagonalise subspace Hamiltonian
!!
!!   (4) Filter by porjecting on trial vectors (localization) 
!!       with weight of Fermi-Dirac functions, f(e) = 1 / ( exp[(e - mu) / kT] + 1 )
!!
!!   This subroutine is called in sub:initial_SFcoeff.
!!
!!  AUTHOR
!!   A.Nakata
!!  CREATION DATE
!!   2017/01/08
!!  MODIFICATION HISTORY
!!   14/11/2019 nakata
!!    Take average of matSFcoeff(1) and matSFcoeff(2) when flag_SpinDependentSF=F and nspin=2
!!
!!  SOURCE
!!
  subroutine LocFilterDiag(LFDhalo)

    use primary_module,         only: bundle
    use cover_module,           only: BCS_parts
    use matrix_module,          only: matrix_halo
    use matrix_data,            only: mat, halo, SFcoeff_range, LD_range
    use mult_module,            only: mult, mat_p, matrix_scale, &
                                      matrix_sum, allocate_temp_matrix, free_temp_matrix, &
                                      matSatomf, matHatomf, matSFcoeff, aLa_aHa_aLHa, aLa_aSa_aLSa
    use global_module,          only: ni_in_cell, numprocs, nspin, nspin_SF, sf, atomf, flag_SpinDependentSF
    use species_module,         only: npao_species
    use io_module,              only: get_file_name
    use input_module,           only: io_assign, io_close
    use GenComms,               only: mtime, gsum
    use timer_stdclocks_module, only: start_timer,stop_timer,tmr_std_allocation

    implicit none

    external :: DSYGVX, dgemm, dcopy
    ! Passed variables
    type(matrix_halo) LFDhalo                         ! LFDhalo = mult(aLa_aHa_aLHa)%ahalo
    ! Local
    integer :: stat, spin, spin_SF
    integer :: lun11,lun12, n_naba_i_d, len_Sub_i_d
    character(len=15) :: filename11,filename12
    real(double) :: ChemP, kT, kT1                   ! filter function
    real(double) :: NEsub, NEsub0                    ! number of electrons in subspace, used for determing ChemP
    integer :: INEsub
    integer :: info                                  ! for DSYGVX
    integer :: matSFcoeff_2
    real(double) :: abstol, t0, t1
    integer :: max_npao, nhalo_LFD, max_npao_LFD, len_kj_sub, len_Sub, &      ! for max subspace of halo atoms
               atom_num, iprim, atom_i, atom_spec, NTVEC,              &      ! for atom_i
               np, i, k, ist, gcspart, k_in_halo, nd3,                 &
               len_Sub_i, LWORK, LIWORK, NUMEIG                               ! for subspace of atom_i 

    integer, allocatable, dimension(:)        :: label_kj_Hsub, label_kj_Ssub ! for subspace of halo atoms
    real(double), allocatable, dimension(:)   :: Hsub, Ssub
    real(double), allocatable, dimension(:,:) :: TVEC, Hsub_i, Ssub_i, EVEC   ! for subspace of atom_i
    real(double), allocatable, dimension(:)   :: EVAL
    integer, allocatable, dimension(:)        :: l_k_g, l_kpao                ! for debug
    real(double), allocatable, dimension(:)   :: l_k_r2
    integer, allocatable, dimension(:)        :: IWORK, IFAIL                 ! scratch space
    real(double), allocatable, dimension(:)   :: WORK, WORK2


    call my_barrier
<<<<<<< HEAD
    if (inode==ionode .and. iprint_basis>2) write(io_lun,*) 'Do Localized Filter Diagonalization'
=======
    if (inode==ionode .and. iprint_basis>1) write(io_lun,*) 'Do Localized Filter Diagonalization'
>>>>>>> Tidy LFD and start LibXC reference output

    do spin_SF = 1, nspin_SF
       call matrix_scale(zero,matSFcoeff(spin_SF))
    enddo

    ! Set variables
    ChemP = LFD_ChemP          ! chemical potential (mu) (default = zero)
    kT    = LFD_kT             ! k*T (default = 0.1)    
    kT1   = one / kT           ! 1 / k*T
    abstol = 1.0e-300_double

    ! open debug file for TVEC and subspace MOs
    if (iprint_basis>=6) then
       call get_file_name('TVECr',numprocs,inode,filename11)  ! Build a filename based on node number for TVEC
       call io_assign(lun11)                                  ! Open file 
       open(unit=lun11,file=filename11)
       call get_file_name('SubMOr',numprocs,inode,filename12) ! Build a filename based on node number for MOs
       call io_assign(lun12)                                  ! Open file
       open(unit=lun12,file=filename12)
    endif

    ! estimate the maximum size of subspace matrices for halo-atoms (Ssub and Hsub)
    max_npao     = maxval(npao_species)                         ! max. number of PAOs belonging to a halo atom
    nhalo_LFD    = LFDhalo%ni_in_halo                           ! number of halo atoms in LD_range
    max_npao_LFD = nhalo_LFD*max_npao                           ! max. number of total PAOs of halo atoms
    len_kj_sub   = (nhalo_LFD*nhalo_LFD+nhalo_LFD)/2            ! max. number of lower-triangle of halo-pair (k,j)
    len_Sub      = (max_npao_LFD*max_npao_LFD+max_npao_LFD)/2   ! max. number of lower-triangle of subspace matrix elements

    if (iprint_basis>=5 .and. inode==ionode) then
       if (flag_LFD_useChemPsub) then 
          write(io_lun,*) ' LFD: Chemical potential = ', 'will be determined later'
       else
          write(io_lun,*) ' LFD: Chemical potential = ', ChemP
       endif
       write(io_lun,*) '      k*T                = ', kT
       write(io_lun,*) '      nhalo_LFD          = ', nhalo_LFD
       write(io_lun,*) '      max_npao_LFD       = ', max_npao_LFD
       write(io_lun,*) '      len_kj_sub         = ', len_kj_sub
       write(io_lun,*) '      len_Sub            = ', len_Sub
    endif

    ! allocate subspace matrices and their halo-pair labels
    call start_timer(tmr_std_allocation)
    allocate(Hsub(len_Sub),Ssub(len_Sub),STAT=stat)
    if(stat/=0) call cq_abort('LocFilterDiag: error allocating H/Ssub')
    call reg_alloc_mem(area_basis, 2*len_Sub, type_dbl)
    allocate(label_kj_Hsub(len_kj_sub),label_kj_Ssub(len_kj_sub),STAT=stat)
    if(stat/=0) call cq_abort('LocFilterDiag: error allocating label_kj_H/Ssub')
    call reg_alloc_mem(area_basis, 2*len_kj_sub, type_int)
    call stop_timer(tmr_std_allocation)

    do spin = 1,nspin
!
!      --- (1) make subspace matrix for halo atoms (NODE)
!
       Hsub(:) = zero
       label_kj_Hsub(:) = 0
       call LFD_make_Subspace_halo(myid,mat_p(matHatomf(spin))%matrix,mat_p(matHatomf(spin))%length,&
                                   Hsub,len_Sub,label_kj_Hsub,len_kj_sub,mult(aLa_aHa_aLHa))
       if (spin.eq.1) then
          Ssub(:) = zero
          label_kj_Ssub(:) = 0
          call LFD_make_Subspace_halo(myid,mat_p(matSatomf)%matrix,mat_p(matSatomf)%length,&
                                      Ssub,len_Sub,label_kj_Ssub,len_kj_sub,mult(aLa_aSa_aLSa))
       endif
!
!      --- Start Localized Filter Diagonalisation Method for each primary ATOM ---
!
       if(myid==0) t0 = mtime()
       if (iprint_basis>=5.and.inode==ionode) &
          write(io_lun,'(/A,I2)') 'Start atom loop in the LFD method for spin',spin

       iprim = 0; atom_i = 0  
       do np = 1,bundle%groups_on_node   ! Loop over primary set partitions
          if(bundle%nm_nodgroup(np)>0) then   ! If there are atoms in partition
             do i = 1,bundle%nm_nodgroup(np)   ! Loop over atom_i
                atom_num = bundle%nm_nodbeg(np) + i - 1
                iprim = iprim + 1
                atom_i = bundle%ig_prim(iprim)   ! global number of i
                atom_spec = bundle%species(atom_num)
                NTVEC = mat(np,SFcoeff_range)%ndimi(i)
!
                ! count the dimension of the subspace for atom_i
                ! can we simplify this?
                len_Sub_i = 0
                do k = 1,mat(np,LD_range)%n_nab(i)                                   ! Loop over atom_k, neighbour of atom_i
                   ist = mat(np,LD_range)%i_acc(i) + k - 1                           ! index of atom_k in neighbour labelling
                   gcspart = BCS_parts%icover_ibeg(mat(np,LD_range)%i_part(ist)) + &
                             mat(np,LD_range)%i_seq(ist) -1                          ! index of atom_k in CS
                   k_in_halo = LFDhalo%i_halo(gcspart)                               ! index of atom_k in halo labelling
                   nd3 = LFDhalo%ndimj(k_in_halo)                                    ! number of PAOs on atom_k
                   len_Sub_i = len_Sub_i + nd3                                       ! dimension of subspace for atom_i
                enddo ! k
                if (iprint_basis>=6) write(io_lun,'(A,I8,A,I2,A,I3,A,I4)') &
                                          '  atom_i=',atom_i, ' : NTVEC=',NTVEC, &
                                          '  nab=',mat(np,LD_range)%n_nab(i), '  len_Sub_i = ',len_Sub_i
!
                ! allocate spaces for atom_i
                call start_timer(tmr_std_allocation)
                allocate(TVEC(len_Sub_i,NTVEC),STAT=stat)                             ! TVEC
                if(stat/=0) call cq_abort('LocFilterDiag: error allocating TVEC')
                allocate(Hsub_i(len_Sub_i,len_Sub_i),STAT=stat)                       ! Hsub_i
                if(stat/=0) call cq_abort('LocFilterDiag: error allocating Hsub_i')
                allocate(Ssub_i(len_Sub_i,len_Sub_i),STAT=stat)                       ! Ssub_i
                if(stat/=0) call cq_abort('LocFilterDiag: error allocating Ssub_i')
                allocate(EVAL(len_Sub_i),STAT=stat)                                   ! EVAL
                if(stat/=0) call cq_abort('LocFilterDiag: error allocating EVAL')
                allocate(EVEC(len_Sub_i,len_Sub_i),STAT=stat)                         ! EVEC
                if(stat/=0) call cq_abort('LocFilterDiag: error allocating EVEC')
                LWORK=max(len_Sub_i*10,len_Sub_i*len_Sub_i)                           ! allocate scratch spaces
                LIWORK=len_Sub_i*5 
                allocate(WORK(LWORK),STAT=stat)                                       ! WORK
                if(stat/=0) call cq_abort('LocFilterDiag: error allocating WORK')
                allocate(WORK2(LWORK),STAT=stat)                                      ! WORK2
                if(stat/=0) call cq_abort('LocFilterDiag: error allocating WORK2')
                allocate(IWORK(LIWORK),STAT=stat)                                     ! IWORK
                if(stat/=0) call cq_abort('LocFilterDiag: error allocating IWORK')
                allocate(IFAIL(len_Sub_i),STAT=stat)                                  ! IFAIL
                if(stat/=0) call cq_abort('LocFilterDiag: error allocating IFAIL')
                ! spaces for saving the information of neighbour atoms for debug
                n_naba_i_d = mat(np,LD_range)%n_nab(i)  
                len_Sub_i_d = len_Sub_i  
                allocate(l_k_g(n_naba_i_d),STAT=stat)                                 ! labels of k for debug
                if(stat/=0) call cq_abort('LocFilterDiag: error allocating l_k_g')
                allocate(l_k_r2(n_naba_i_d),STAT=stat)                                ! labels of k for debug
                if(stat/=0) call cq_abort('LocFilterDiag: error allocating l_k_r2')
                allocate(l_kpao(len_Sub_i_d),STAT=stat)                               ! labels of k_pao for debug
                if(stat/=0) call cq_abort('LocFilterDiag: error allocating l_kpao')
                call stop_timer(tmr_std_allocation)
                ! zero clear
                TVEC(:,:)    = zero
                Hsub_i(:,:)  = zero
                Ssub_i(:,:)  = zero
                EVAL(:)      = zero
                EVEC(:,:)    = zero
                WORK(:)      = zero
                WORK2(:)     = zero
                IWORK(:)     = 0
                IFAIL(:)     = 0
                l_k_g(:)     = 0
                l_k_r2(:)    = zero
                l_kpao(:)    = 0
!
!               --- make trial vectors for atom_i ---
!
                NEsub = zero
                call LFD_make_TVEC(TVEC,NTVEC,len_Sub_i,np,i,atom_i,atom_num,atom_spec,LFDhalo, NEsub, &
                                   len_Sub_i_d,n_naba_i_d,l_k_g,l_k_r2,l_kpao)  
                !call LFD_debug_matrix(lun11,0,np,atom_i,EVAL,TVEC,len_Sub_i,NTVEC, &
                !                                           n_naba_i_d,l_k_g,l_k_r2,l_kpao)
!
!               --- (2) make subspaces for atom_i
!
                call LFD_make_Subspace_i(np,i,LFDhalo,Ssub,Hsub,label_kj_Ssub,label_kj_Hsub, &
                                         Ssub_i,Hsub_i,len_Sub_i,len_Sub,len_kj_sub)
!                write(io_lun,'(/A,I7)') '  Hsub_i:',atom_i
!                call LFD_debug_matrix(io_lun,0,np,atom_i,EVAL,Hsub_i,len_Sub_i,len_Sub_i, &
!                                      n_naba_i_d,l_k_g,l_k_r2,l_kpao)
!                write(io_lun,'(/A,I7)') '  Ssub_i:',atom_i
!                call LFD_debug_matrix(io_lun,0,np,atom_i,EVAL,Ssub_i,len_Sub_i,len_Sub_i, &
!                                      n_naba_i_d,l_k_g,l_k_r2,l_kpao)
!
!               --- (3) diagonalise subspace matrix: Hsub_i * Csub_i = e * Ssub_i * Csub_i
!
                call DCOPY(len_Sub_i*len_Sub_i,Ssub_i,1,WORK2,1)     ! save original Ssub_i
                call DSYGVX(1, 'V', 'A', 'U', len_Sub_i, Hsub_i, len_Sub_i, Ssub_i, len_Sub_i, &
                            0.0d0, 0.0d0, 0, 0, abstol, NUMEIG, EVAL, EVEC, len_Sub_i, WORK,   &
                            LWORK, IWORK, IFAIL, INFO)
                if (NUMEIG.ne.len_Sub_i) &
                   call cq_abort('In sub: LocFilterDiag, error in diagonalization for atom ',atom_i)
                Ssub_i(:,:) = zero
                call DCOPY(len_Sub_i*len_Sub_i,WORK2,1,Ssub_i,1)     ! Ssub_i was overwritten, so reset to the original.
                ! debug: print out local MOs
                if (iprint_basis>=6) call LFD_debug_matrix(lun12,1,np,atom_i,EVAL,EVEC,len_Sub_i,NUMEIG, &
                                                           n_naba_i_d,l_k_g,l_k_r2,l_kpao)
!
!               --- (4) filteration: k = Csub_i * f(e) * Csub_i**T * Ssub_i * TVEC
!
!               (4-1) t~ = Csub_i**TVEC * Ssub_i * t, on exit WORK contains t~ and copied into TVEC
                WORK(:) = zero
                call DGEMM('N','N',len_Sub_i,NTVEC,len_Sub_i,one,Ssub_i,len_Sub_i,TVEC,len_Sub_i,zero,WORK,len_Sub_i)
                TVEC(:,:) = zero
                call DCOPY(len_Sub_i*NTVEC,WORK,1,TVEC,1)
!
                WORK(:) = zero
                call DGEMM('T','N',NUMEIG,NTVEC,len_Sub_i,one,EVEC,len_Sub_i,TVEC,len_Sub_i,zero,WORK,NUMEIG)
                TVEC(:,:) = zero
                call DCOPY(NUMEIG*NTVEC,WORK,1,TVEC,1)
!
!               (4-2) multiply filtration function f(e): t~ = f(e) * t~, on exit, TVEC contains t~.
                if (flag_LFD_useChemPsub) then
                   NEsub0 = NEsub / two                          ! NEsub is the number of electrons in the subspace
                   NEsub  = dint(NEsub0) 
                   if (NEsub0.gt.NEsub) NEsub = NEsub + one
                   INEsub = idint(NEsub)                         ! number of occupied orbitals
                   ChemP = (EVAL(INEsub)+EVAL(INEsub+1)) / two   ! ChemP is set to the average of subspace HOMO and LUMO
                   if (iprint_basis>=6) write(io_lun,'(A,i8,A,F10.5)') ' Atom',atom_i,': ChemP of this subspace =',ChemP
                endif
                call LFD_filter(TVEC,EVAL,ChemP,kT1,NUMEIG,NTVEC)
!
!               (4-3) t~ = d * t~, on exit WRK contains t~ and copied into TVEC
                WORK(:) = zero
                call DGEMM('N','N',len_Sub_i,NTVEC,NUMEIG,one,EVEC,len_Sub_i,TVEC,NUMEIG,zero,WORK,len_Sub_i)
                TVEC(:,:) = zero
                call DCOPY(len_Sub_i*NTVEC,WORK,1,TVEC,1)
!          
!               --- copy filtered vectors into the corresponding positions in matSFcoeff ---
!
                WORK(:) = zero
                call transpose_2Dmat(TVEC,WORK,len_Sub_i,NTVEC)
                if (.not.flag_SpinDependentSF .and. spin.eq.2) then
                   if(inode==ionode  .and. iprint_basis>1) &
                        write(io_lun,*) 'Take average of matSFcoeff(1) and matSFcoeff(2) into matSFcoeff(1)'
                   matSFcoeff_2 = allocate_temp_matrix(SFcoeff_range,0,sf,atomf)
                   call matrix_scale(zero,matSFcoeff_2)
                   call LFD_put_TVEC_to_SFcoeff(np,i,mat(:,SFcoeff_range),mat_p(matSFcoeff_2)%matrix, &
                                                mat_p(matSFcoeff_2)%length,mat(:,LD_range),WORK,len_Sub_i*NTVEC)
                   call matrix_sum(half, matSFcoeff(1), half, matSFcoeff_2)
                   call free_temp_matrix(matSFcoeff_2)
                else
                   call LFD_put_TVEC_to_SFcoeff(np,i,mat(:,SFcoeff_range),mat_p(matSFcoeff(spin))%matrix, &
                                                mat_p(matSFcoeff(spin))%length,mat(:,LD_range),WORK,len_Sub_i*NTVEC)
                endif
!
!               --- deallocate spaces for atom_i ---
!                                  
                call start_timer(tmr_std_allocation)                                  
                deallocate(l_kpao,STAT=stat)
                if(stat/=0) call cq_abort('LocFilterDiag: error deallocating l_kpao')
                deallocate(l_k_r2,STAT=stat)
                if(stat/=0) call cq_abort('LocFilterDiag: error deallocating l_k_r2')
                deallocate(l_k_g,STAT=stat)
                if(stat/=0) call cq_abort('LocFilterDiag: error deallocating l_k_g')
                deallocate(IFAIL,STAT=stat)
                if(stat/=0) call cq_abort('LocFilterDiag: error deallocating IFAIL')
                deallocate(IWORK,STAT=stat)
                if(stat/=0) call cq_abort('LocFilterDiag: error deallocating IWORK')
                deallocate(WORK2,STAT=stat)
                if(stat/=0) call cq_abort('LocFilterDiag: error deallocating WORK2')
                deallocate(WORK,STAT=stat)
                if(stat/=0) call cq_abort('LocFilterDiag: error deallocating WORK')
                deallocate(EVEC,STAT=stat)
                if(stat/=0) call cq_abort('LocFilterDiag: error deallocating EVEC')
                deallocate(EVAL,STAT=stat)
                if(stat/=0) call cq_abort('LocFilterDiag: error deallocating EVAL')
                deallocate(Ssub_i,STAT=stat)
                if(stat/=0) call cq_abort('LocFilterDiag: error deallocating Ssub_i')
                deallocate(Hsub_i,STAT=stat)
                if(stat/=0) call cq_abort('LocFilterDiag: error deallocating Hsub_i')
                deallocate(TVEC,STAT=stat)
                if(stat/=0) call cq_abort('LocFilterDiag: error deallocating TVEC')
                call stop_timer(tmr_std_allocation)
!             
             enddo ! i
          endif ! endif of (bundle%nm_nodgroup(np)>0)
       enddo ! np
    enddo ! spin

    ! deallocate subspace matrices their labels
    call start_timer(tmr_std_allocation)
    deallocate(label_kj_Hsub,label_kj_Ssub,STAT=stat)
    if(stat/=0) call cq_abort('LocFilterDiag: error deallocating label_kj_H/Ssub')
    call reg_dealloc_mem(area_basis, 2*len_kj_sub, type_int)
    deallocate(Hsub,Ssub,STAT=stat)
    if(stat/=0) call cq_abort('LocFilterDiag: error deallocating H/Ssub')
    call reg_dealloc_mem(area_basis, 2*len_Sub, type_dbl)
    call stop_timer(tmr_std_allocation)

    ! Close debug file
    if (iprint_basis>=6) then
       call io_close(lun11)
       call io_close(lun12)
    endif

<<<<<<< HEAD
    if(myid==0 .and. iprint_basis>3) then
       t1 = mtime()
       write(io_lun,'(A,f20.8,A)') 'Time for Localised Filter Diagonalisation: ',t1-t0,' ms'
    end if
    if (inode==ionode .and. iprint_basis>3) write(io_lun,*) 'Done Localized Filter Diagonalization'
=======
    if(myid==0 .and. inode==ionode  .and. iprint_basis>1) then
       t1 = mtime()
       write(io_lun,'(A,f20.8,A)') 'Time for Localised Filter Diagonalisation: ',t1-t0,' ms'
    end if
    if (inode==ionode .and. iprint_basis>2) write(io_lun,*) 'Done Localized Filter Diagonalization'
>>>>>>> Tidy LFD and start LibXC reference output
!
    return
  end subroutine LocFilterDiag
!!***

!!****f* multisiteSF_module/LFD_make_Subspace_halo *
!!
!!  NAME 
!!   LFD_make_Subspace_halo
!! 
!!  PURPOSE
!!   This subroutine makes the subspaces of halo atoms for this NODE.
!!
!!   The size of the subspace matrix for atom_i is
!!        (neighbour of atom_i)*(neighbour of atom_i).
!!   Therefore, the size of the subspace for this NODE is
!!        (halo atoms for LD_range)*(halo atoms for LD_range).
!!
!!   Get the information of matrix B from remote nodes
!!   and make subspace Bsub(halo,halo) in lower-triangle form
!!   by calling sub:LFD_pickup_Subspace_elements.
!!
!!   Bsub contains only lower-triangle elements of the subspace matrix
!!   and only the elements within Brange.
!!   The label of the elements is saved in label_kj.
!!
!!   This subroutine is based on sub:mat_mult.
!!
!!   This subroutine is called in sub:LocFilterDiag.
!!
!!  AUTHOR
!!   A.Nakata
!!  CREATION DATE
!!   2017/01/08
!!  MODIFICATION HISTORY
!!   2018/10/08 07:59 dave
!!    Adding recv_part variable to count partitions received from processes
!!    to follow update of multiplies (to conform to MPI standard for
!!    tags)
!!
!!  SOURCE
!!
  subroutine LFD_make_Subspace_halo(myid,b,lenb,bsub,len_bsub,label_kj_sub,len_kj_sub,a_b_c)

    use global_module
    use matrix_module
    use multiply_module, ONLY: prefetch
    use basic_types
    use matrix_comms_module
    use comms_module
    use GenComms, ONLY: mtime
    use mpi

    implicit none

    ! Mult type - contains indices, small groups etc
    type(matrix_mult)::a_b_c
    ! Matrices
    integer :: lenb, len_bsub, len_kj_sub
    real(double) :: b(lenb)
    real(double) :: bsub(len_bsub)
    integer :: label_kj_sub(len_kj_sub)

    ! Node id - should be in a module ?
    integer :: myid

    ! Local variables
    ! This will be dynamically allocated/deallocated by the system
    integer :: lab_const
    integer :: ierr,kpart,ind_part,ncover_yz,n_which,ipart,nnode
    integer :: icall,n_cont,kpart_next,ind_partN,k_off
    integer :: stat,ilen2,lenb_rem
    ! Remote variables to be allocated
    integer(integ),allocatable :: ibpart_rem(:)
    real(double),allocatable :: b_rem(:)
    ! Remote variables which will point to part_array
    integer(integ),pointer :: nbnab_rem(:)
    integer(integ),pointer :: ibseq_rem(:)
    integer(integ),pointer :: ibind_rem(:)
    integer(integ),pointer :: ib_nd_acc_rem(:)
    integer(integ),pointer :: npxyz_rem(:)
    integer(integ),pointer :: ibndimj_rem(:)
    ! Arrays for remote variables to point to
    integer, target :: part_array(3*a_b_c%parts%mx_mem_grp+ &
         5*a_b_c%parts%mx_mem_grp*a_b_c%bmat(1)%mx_abs)
    integer, dimension(:), allocatable :: nreqs
    integer :: offset,sends,i,j
    integer, dimension(MPI_STATUS_SIZE) :: mpi_stat

    integer, allocatable, dimension(:) :: recv_part

    logical flag,call_flag
    real(double) :: t0,t1

    integer :: count_kj


    if (iprint_basis>=5 .and. inode==ionode) write(io_lun,*) 'We are in sub:LFD_make_Subspace_halo'

    call start_timer(tmr_std_allocation)
    if(iprint_mat>3.AND.myid==0) t0 = mtime()    
    ! Allocate memory for the elements
    allocate(ibpart_rem(a_b_c%parts%mx_mem_grp*a_b_c%bmat(1)%mx_abs),STAT=stat)
    if(stat/=0) call cq_abort('LFD_make_Subspace_halo: error allocating ibpart_rem')
    allocate(recv_part(0:a_b_c%comms%inode),STAT=stat)
    if(stat/=0) call cq_abort('LFD_make_Subspace_halo: error allocating recv_part')
    recv_part = zero
    call stop_timer(tmr_std_allocation)

    sends = 0
    do i=1,a_b_c%comms%inode
       if(a_b_c%comms%np_send(i)>0) then
          do j=1,a_b_c%comms%np_send(i)
             sends = sends+1
          end do
       end if
    end do

    call start_timer(tmr_std_allocation)
    allocate(nreqs(sends*2),STAT=stat)
    if(stat/=0) call cq_abort("LFD_make_Subspace_halo: Error allocating nreqs",sends,stat)
    call stop_timer(tmr_std_allocation)

    sends = 0
    ! zero subspace matrix and its label
    label_kj_sub(:) = 0
    Bsub(:) = zero
    count_kj = 0

    ! Start the send procedure
    !write(io_lun,*) 'Sending'
    call Mquest_start_send(a_b_c,b,nreqs,myid,a_b_c%prim%mx_ngonn,sends)
    !write(io_lun,*) 'Returned ',a_b_c%ahalo%np_in_halo,myid
    ncover_yz=a_b_c%gcs%ncovery*a_b_c%gcs%ncoverz

    do kpart = 1,a_b_c%ahalo%np_in_halo  ! Main loop. Loop for halo-partition kpart
       !write(io_lun,*) 'Part: ',kpart,myid
       icall=1
       ind_part = a_b_c%ahalo%lab_hcell(kpart)
       !write(io_lun,*) 'ind_part: ',ind_part

       if(kpart>1) then  ! Is it a periodic image of the previous partition ?
          if(ind_part.eq.a_b_c%ahalo%lab_hcell(kpart-1)) then 
             icall=0
          else ! Get the data
             !write(io_lun,*) myid,' seq: ',size(a_b_c%parts%i_cc2seq)
             ipart = a_b_c%parts%i_cc2seq(ind_part)
             !write(io_lun,*) myid,' Alloc b_rem part: ',ipart
             nnode = a_b_c%comms%neigh_node_list(kpart)
             recv_part(nnode) = recv_part(nnode)+1
             !write(io_lun,*) myid,' Alloc b_rem node: ',nnode
             !write(io_lun,*) myid,' Alloc b_rem icc: ',a_b_c%parts%i_cc2node(ind_part)
             !write(io_lun,*) myid,' Alloc b_rem alloc: ',allocated(b_rem)
             if(allocated(b_rem)) deallocate(b_rem)
             if(a_b_c%parts%i_cc2node(ind_part)==myid+1) then
                lenb_rem = a_b_c%bmat(ipart)%part_nd_nabs
             else
                lenb_rem = a_b_c%comms%ilen3rec(ipart,nnode)
             end if
             allocate(b_rem(lenb_rem))
             call prefetch(kpart,a_b_c%ahalo,a_b_c%comms,a_b_c%bmat,icall,&  
                  n_cont,part_array,a_b_c%bindex,b_rem,lenb_rem,b,myid,ilen2,&
                  mx_msg_per_part,a_b_c%parts,a_b_c%prim,a_b_c%gcs,(recv_part(nnode)-1)*2)
             !write(io_lun,*) 'b_rem: ',lenb_rem
             ! Now point the _rem variables at the appropriate parts of 
             ! the array where we will receive the data
             offset = 0
             nbnab_rem => part_array(offset+1:offset+n_cont)
             offset = offset+n_cont
             ibind_rem => part_array(offset+1:offset+n_cont)
             offset = offset+n_cont
             ib_nd_acc_rem => part_array(offset+1:offset+n_cont)
             offset = offset+n_cont
             ibseq_rem => part_array(offset+1:offset+ilen2)
             offset = offset+ilen2
             npxyz_rem => part_array(offset+1:offset+3*ilen2)
             offset = offset+3*ilen2
             ibndimj_rem => part_array(offset+1:offset+ilen2)
             if(offset+ilen2>3*a_b_c%parts%mx_mem_grp+ &
                  5*a_b_c%parts%mx_mem_grp*a_b_c%bmat(1)%mx_abs) then
                call cq_abort('mat_mult: error pointing to part_array ',kpart)
             endif
             ! Create ibpart_rem
             call end_part_comms(myid,n_cont,nbnab_rem,ibind_rem,npxyz_rem,&
                  ibpart_rem,ncover_yz,a_b_c%gcs%ncoverz)
          endif
       else ! Get the data
          !write(io_lun,*) myid,' seq: ',size(a_b_c%parts%i_cc2seq)
          ipart = a_b_c%parts%i_cc2seq(ind_part)
          !write(io_lun,*) myid,' Alloc b_rem part: ',ipart
          nnode = a_b_c%comms%neigh_node_list(kpart)
          recv_part(nnode) = recv_part(nnode)+1
          !write(io_lun,*) myid,' Alloc b_rem node: ',nnode
          !write(io_lun,*) myid,' Alloc b_rem icc: ',a_b_c%parts%i_cc2node(ind_part)
          !write(io_lun,*) myid,' Alloc b_rem alloc: ',allocated(b_rem)
          if(allocated(b_rem)) deallocate(b_rem)
          if(a_b_c%parts%i_cc2node(ind_part)==myid+1) then
             lenb_rem = a_b_c%bmat(ipart)%part_nd_nabs
          else
             lenb_rem = a_b_c%comms%ilen3rec(ipart,nnode)
          end if
          call start_timer(tmr_std_allocation)
          allocate(b_rem(lenb_rem))
          call stop_timer(tmr_std_allocation)
          call prefetch(kpart,a_b_c%ahalo,a_b_c%comms,a_b_c%bmat,icall,& 
               n_cont,part_array,a_b_c%bindex,b_rem,lenb_rem,b,myid,ilen2,&
               mx_msg_per_part,a_b_c%parts,a_b_c%prim,a_b_c%gcs,(recv_part(nnode)-1)*2)
          lenb_rem = size(b_rem)
          !write(io_lun,*) 'b_rem: ',lenb_rem
          ! Now point the _rem variables at the appropriate parts of the array 
          ! where we will receive the data
          offset = 0
          nbnab_rem => part_array(offset+1:offset+n_cont)
          offset = offset+n_cont
          ibind_rem => part_array(offset+1:offset+n_cont)
          offset = offset+n_cont
          ib_nd_acc_rem => part_array(offset+1:offset+n_cont)
          offset = offset+n_cont
          ibseq_rem => part_array(offset+1:offset+ilen2)
          offset = offset+ilen2
          npxyz_rem => part_array(offset+1:offset+3*ilen2)
          offset = offset+3*ilen2
          ibndimj_rem => part_array(offset+1:offset+ilen2)
          if(offset+ilen2>3*a_b_c%parts%mx_mem_grp+ &
               5*a_b_c%parts%mx_mem_grp*a_b_c%bmat(1)%mx_abs) then
             call cq_abort('Error pointing to part_array !',kpart)
          endif
          call end_part_comms(myid,n_cont,nbnab_rem,ibind_rem,npxyz_rem,&
               ibpart_rem,ncover_yz,a_b_c%gcs%ncoverz)
       endif ! End of the "if this isn't the first partition" loop
       k_off=a_b_c%ahalo%lab_hcover(kpart) ! --- offset for pbcs

       call LFD_pickup_Subspace_elements(k_off,kpart,ib_nd_acc_rem,ibind_rem,nbnab_rem,&
                                         ibpart_rem,ibseq_rem,ibndimj_rem,&
                                         b_rem,a_b_c%ahalo,a_b_c%chalo,a_b_c%ltrans,&
                                         a_b_c%bmat(1)%mx_abs,a_b_c%parts%mx_mem_grp, &
                                         lenb_rem,label_kj_sub,Bsub,len_kj_sub,len_bsub,count_kj)

    enddo ! End of the kpart=1,ahalo%np_in_halo loop !

    call start_timer(tmr_std_allocation)
    if(allocated(b_rem)) deallocate(b_rem)
    call stop_timer(tmr_std_allocation)

    ! --------------------------------------------------
    ! End of the main loop over partitions in the A-halo
    ! --------------------------------------------------
    ! sends is only incremented in the MPI version of Mquest_start_send, so this works
    !write(io_lun,*) 'Done loop ',myid,sends
    if(sends>0) then
       do i=1,sends
          call MPI_Wait(nreqs(i),mpi_stat,ierr)
          if(ierr/=0) call cq_abort("Error waiting for send to finish",i)
          !write(io_lun,*) 'Send done ',i,myid
       end do
    end if
    call my_barrier
    call start_timer(tmr_std_allocation)
    deallocate(nreqs,STAT=stat)
    if(stat/=0) call cq_abort('LFD_make_Subspace_halo: error deallocating nreqs')
    deallocate(ibpart_rem,STAT=stat)
    if(stat/=0) call cq_abort('LFD_make_Subspace_halo: error deallocating ibpart_rem')
    deallocate(recv_part,STAT=stat)
    if(stat/=0) call cq_abort('LFD_make_Subspace_halo: error deallocating recv_part')
    !deallocate(b_rem,STAT=stat)
    !if(stat/=0) call cq_abort('LFD_make_Subspace_halo: error deallocating b_rem')
    call stop_timer(tmr_std_allocation)
    call my_barrier
    if(iprint_mat>3.AND.myid==0) then
       t1 = mtime()
       write(io_lun,'(A,f20.8,A)') 'Time for LFD_make_Subspace_halo: ',t1-t0,' ms'
    end if
!
    return
  end subroutine LFD_make_Subspace_halo
!!***

!!****f* multisiteSF_module/LFD_pickup_Subspace_elements *
!!
!!  NAME 
!!   LFD_pickup_Subspace_elements
!! 
!!  PURPOSE
!!   Matrix elements of B (= Hpao or Spao) belonging to neighbour atoms 
!!   within LD_range will be put into subspace matrix Bsub.
!!   The indices (k,j) and the matrix elements belonging to these atoms 
!!   are saved into label_kj_sub and Bsub.
!!
!!
!!   k_off:          offset in partition labelling to account for p.b.c.
!!
!!   kpart:          A-halo seq. no. of current partition K
!!
!!   ib_nd_acc(k_in_part): address in arrays ibseq(.) where data for atom_k start.
!!
!!   ibaddr(k_in_part):    address in arrays ibpart(.) where matrix element data for atom_k start.
!!
!!   nbnab(k_in_part):     number of B-neighbours of atom k, the latter being
!!                         specified by its (unpruned) seq. no. in partition K.
!!
!!   ibpart(.):      partition to which each atom j belongs, the latter
!!                   being specified as jth B-neighbour of atom k.
!!
!!   ibseq(.):       unpruned sequence number in partition for each
!!                   atom j, the latter being specified as in ibpart.
!!
!!   bndim2(.):      number of the 2nd dimension of B belonging to atom_j
!!
!!   This subroutine is based on sub:m_kern_max.
!!
!!   This subroutine is called in sub:LFD_make_Subspace_halo.
!!
!!  AUTHOR
!!   A.Nakata
!!  CREATION DATE
!!   2017/01/08
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine LFD_pickup_Subspace_elements(k_off,kpart,&
       ib_nd_acc,ibaddr,nbnab,ibpart,ibseq,bndim2,b,ahalo,chalo,at,&
       mx_absb,mx_part,lenb,label_kj_sub,Bsub,len_kj_sub,len_bsub,count_kj)

    use matrix_module
    use basic_types,    ONLY: primary_set
    use primary_module, ONLY: bundle

    implicit none

    ! Passed variables
    type(matrix_halo)::ahalo,chalo
    type(matrix_trans)::at
    integer :: mx_absb,mx_part,lenb
    real(double) :: b(lenb)
    integer :: len_kj_sub,len_bsub
    integer :: label_kj_sub(len_kj_sub)
    real(double) :: Bsub(len_bsub)
    integer :: count_kj

    integer :: kpart,k_off
    ! Remote indices
    integer(integ) :: ib_nd_acc(mx_part)
    integer(integ) :: ibaddr(mx_part)
    integer(integ) :: nbnab(mx_part)
    integer(integ) :: ibpart(mx_part*mx_absb)
    integer(integ) :: ibseq(mx_part*mx_absb)
    integer(integ) :: bndim2(mx_part*mx_absb)
       
    ! Local variables
    integer :: jbnab2ch(mx_absb)  ! Automatic array
    integer :: nbkbeg,k,k_in_part,k_in_halo,j,jpart,jseq
    integer :: i,i_in_prim,icad,nbbeg,j_in_halo,ncbeg
    integer :: n2,n3,nb_nd_kbeg
    integer :: nd2,nd3
    integer :: nbaddr
!    integer :: kj, count_kj
    integer :: kj


!    if (inode==ionode) write(io_lun,*) 'We are in sub:LFD_pickup_Subspace_elements'

    kj = 0

    ! Loop over halo-atom k in current A-halo partn (=kpart)
    do k=1,ahalo%nh_part(kpart)
       k_in_halo=ahalo%j_beg(kpart)+k-1       ! halo-index of atom_k for A with LD_range
       k_in_part=ahalo%j_seq(k_in_halo)       ! index of atom_k in kpart
       nbkbeg=ibaddr(k_in_part)               ! where starts the information of the neighbour atoms of atom_k for B
       nb_nd_kbeg=ib_nd_acc(k_in_part)        ! where starts the matrix elements of atom_k for B
       nd3 = ahalo%ndimj(k_in_halo)           ! number of PAOs belonging to atom_k

       ! transcription of j from partition to C-halo labelling
       ! (j = neighbour atoms of k for B with Brange)
       do j=1,nbnab(k_in_part)
          jpart=ibpart(nbkbeg+j-1)+k_off
          jseq=ibseq(nbkbeg+j-1)
          jbnab2ch(j)=chalo%i_halo(chalo%i_hbeg(jpart)+jseq-1)   ! index of atom_j in C-halo labelling 
       enddo

       ! Loop over atom_i (atoms in primary set, of which atom_k is a A-neighbour with LD_range)
       do i=1,at%n_hnab(k_in_halo)
          i_in_prim=at%i_prim(at%i_beg(k_in_halo)+i-1)
          icad=(i_in_prim-1)*chalo%ni_in_halo
          nbbeg=nb_nd_kbeg   ! starting position of the B-matrix elements of atom_k

          ! Loop over atom_j (B-neighbour of atom_k)
          do j=1,nbnab(k_in_part)
             nd2 = bndim2(nbkbeg+j-1)   ! number of PAOs on atom_j
             j_in_halo=jbnab2ch(j)      ! index of atom_j in C-halo labelling

             if(j_in_halo.ne.0) then    ! if atom_j is a halo-atom of atom_i for C with LD_range
                ncbeg=chalo%i_h2d(icad+j_in_halo)
                if(ncbeg/=0) then       ! if atom_j is a neighbour atom of atom_i for C with LD_range

                   ! neighbour-atom loops should be canonical (i.e., only j =< k is saved)
                   if (j_in_halo .lt. k_in_halo) then
                      kj = (k_in_halo*k_in_halo - k_in_halo)/2 + j_in_halo   ! canonical number of k and j
                      if (label_kj_sub(kj).eq.0) then   ! if the combination (k,j) didn't appear yet
                         label_kj_sub(kj) = count_kj + 1   ! starting position of PAO-PAO matrix elements for atoms (k,j)
                         do n2=1,nd2   ! loop over PAOs on atom_j (= column of Hpao)
                            nbaddr = nbbeg+nd3*(n2-1)
                            do n3=1,nd3   ! loop over PAOs on atom_k (= row of Hpao)
                               count_kj = count_kj + 1
                               Bsub(count_kj) = b(nbaddr+n3-1)   ! b = Hpao or Spao
                            enddo ! n3
                         enddo ! n2
                      endif ! End of if (label_kj_sub.ne.0)
                   else if (j_in_halo .eq. k_in_halo) then   ! if k = j, PAO loops also should be canonical.
                      kj = (k_in_halo*k_in_halo - k_in_halo)/2 + j_in_halo   ! canonical number of k and j
                      if (label_kj_sub(kj).eq.0) then   ! if the combination (k,j) didn't appear yet
                         label_kj_sub(kj) = count_kj + 1   ! starting position of PAO-PAO matrix elements for atoms (k,j)
                         do n2=1,nd2   ! loop over PAOs on atom_j (= column of Hpao)
                            nbaddr = nbbeg+nd3*(n2-1)
                            do n3=1,n2   ! loop over PAOs on atom_k (= row of Hpao)
                               count_kj = count_kj + 1
                               Bsub(count_kj) = b(nbaddr+n3-1)   ! b = Hpao or Spao
                            enddo ! n3
                         enddo ! n2
                      endif ! End of if (label_kj_sub.ne.0)
                   endif ! End of if (j_in_halo.le.k_in_halo)
                   
                endif ! End of if (ncbeg.ne.0)
             endif ! End of if (j_in_halo.ne.0)

             nbbeg = nbbeg + nd3*nd2
          enddo ! j
       enddo ! i
    enddo ! k
!    
    return
  end subroutine LFD_pickup_Subspace_elements
!!***

!!****f* multisiteSF_module/LFD_make_Subspace_i *
!!
!!  NAME 
!!   LFD_make_Subspace_i
!! 
!!  PURPOSE
!!   This subroutine makes the subspaces for this ATOM.
!!
!!   The size of the subspace matrix for atom_i is
!!       (PAOs of neighbours of atom_i)*(PAOs of neighbours of atom_i).
!!
!!   The subspaces for this NODE were made already by sub:LFD_make_Subspace_halo
!!   and were saved in Ssub and Hsub.
!!   Therefore, the elements for atom_i will be picked up 
!!   and be put into Ssub_i and Hsub_i in this subroutine.
!!
!!   This subroutine is based on sub:m_kern_max and sub:assemble_2.
!!
!!   This subroutine is called in sub:LocFilterDiag.
!!
!!  AUTHOR
!!   A.Nakata
!!  CREATION DATE
!!   2017/01/08
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine LFD_make_Subspace_i(np,i,LFDhalo,Ssub,Hsub,label_kj_Ssub,label_kj_Hsub, &
                                 Ssub_i,Hsub_i,len_Sub_i,len_Sub,len_kj_sub)
                            
    use cover_module,  ONLY: BCS_parts
    use matrix_data,   ONLY: mat, LD_range
    use matrix_module, ONLY: matrix_halo

    implicit none

    ! Passed variables
    type(matrix_halo)  :: LFDhalo
    integer, intent(in) :: np,i,len_Sub_i,len_Sub,len_kj_sub
    real(double), intent(in), dimension(len_Sub) :: Ssub,Hsub
    integer, intent(in), dimension(len_kj_sub) :: label_kj_Ssub,label_kj_Hsub
    real(double), dimension(len_Sub_i,len_Sub_i) :: Ssub_i,Hsub_i

    ! Local variables
    integer :: k, kpao, kpao0, jpao, jpao0, kj, naddr, &
               ist, gcspart, k_in_halo, j, ist2, gcspart2, j_in_halo, nd2, nd3, n2, n3


<<<<<<< HEAD
    if (inode==ionode .and. iprint_basis>=5) write(io_lun,*) 'We are in sub:LFD_make_Subspace_i'
=======
    if (inode==ionode .and. iprint_basis>1) write(io_lun,*) 'We are in sub:LFD_make_Subspace_i'
>>>>>>> Tidy LFD and start LibXC reference output

    kpao0 = 0
    ! Loop over atom_k, neighbour of atom_i
    do k = 1,mat(np,LD_range)%n_nab(i)
       ist = mat(np,LD_range)%i_acc(i) + k - 1                           ! index of atom_k in neighbour labelling
       gcspart = BCS_parts%icover_ibeg(mat(np,LD_range)%i_part(ist)) + &
                 mat(np,LD_range)%i_seq(ist) -1                          ! index of atom_k in CS
       k_in_halo = LFDhalo%i_halo(gcspart)                               ! index of atom_k in halo labelling
       nd3 = LFDhalo%ndimj(k_in_halo)

       jpao0 = 0
       ! Loop over atom_j, neighbour of atom_i
       do j = 1,mat(np,LD_range)%n_nab(i)
          ist2 = mat(np,LD_range)%i_acc(i) + j - 1                            ! index of atom_j in neighbour labelling
          gcspart2 = BCS_parts%icover_ibeg(mat(np,LD_range)%i_part(ist2)) + &
                     mat(np,LD_range)%i_seq(ist2) -1                          ! index of atom_j in CS
          j_in_halo = LFDhalo%i_halo(gcspart2)                                ! index of atom_j in halo labelling
          nd2 = LFDhalo%ndimj(j_in_halo)
  
          ! neighbour-atom loops should be lower-triangle (i.e., only j =< k is saved)
          if (j_in_halo .lt. k_in_halo) then
             kj = (k_in_halo*k_in_halo - k_in_halo)/2 + j_in_halo   ! canonical index of the pair (halo-k, halo-j)
             if (label_kj_Hsub(kj).ne.0) then   ! matrix elements for (k,j) were saved in Hsub
                naddr = label_kj_Hsub(kj) - 1
                jpao = jpao0
                do n2=1,nd2   ! loop over PAOs on atom_j (= column of Hpao)
                   jpao = jpao + 1
                   kpao = kpao0
                   do n3=1,nd3   ! loop over PAOs on atom_k (= row of Hpao)
                      naddr = naddr + 1
                      kpao = kpao + 1
                      Hsub_i(kpao,jpao) = Hsub(naddr)
                   enddo ! n3
                enddo ! n2
             endif ! End of if (label_kj_Hsub.ne.0)
             if (label_kj_Ssub(kj).ne.0) then   ! matrix elements for (k,j) were saved in Ssub
                naddr = label_kj_Ssub(kj) - 1
                jpao = jpao0
                do n2=1,nd2   ! loop over PAOs on atom_j (= column of Spao)
                   jpao = jpao + 1
                   kpao = kpao0
                   do n3=1,nd3   ! loop over PAOs on atom_k (= row of Spao)
                      naddr = naddr + 1
                      kpao = kpao + 1
                      Ssub_i(kpao,jpao) = Ssub(naddr)
                   enddo ! n3
                enddo ! n2
             endif ! End of if (label_kj_Ssub.ne.0)
          else if (j_in_halo .eq. k_in_halo) then   ! if k = j, PAO loops also should be canonical.
             kj = (k_in_halo*k_in_halo - k_in_halo)/2 + j_in_halo   ! canonical index of the pair (halo-k,halo-j)
             if (label_kj_Hsub(kj).ne.0) then   ! matrix elements for (k,j) were saved in Hsub
                naddr = label_kj_Hsub(kj) - 1
                jpao = jpao0
                do n2=1,nd2   ! loop over PAOs on atom_j (= column of Hpao)
                   jpao = jpao + 1
                   kpao = kpao0
                   do n3=1,n2   ! loop over PAOs on atom_k (= row of Hpao)
                      naddr = naddr + 1
                      kpao = kpao + 1
                      Hsub_i(kpao,jpao) = Hsub(naddr)
                   enddo ! n3
                enddo ! n2
             endif ! End of if (label_kj_Hsub.ne.0)
             if (label_kj_Ssub(kj).ne.0) then   ! matrix elements for (k,j) were saved in Ssub
                naddr = label_kj_Ssub(kj) - 1
                jpao = jpao0
                do n2=1,nd2   ! loop over PAOs on atom_j (= column of Spao)
                   jpao = jpao + 1
                   kpao = kpao0
                   do n3=1,n2   ! loop over PAOs on atom_k (= row of Spao)
                      naddr = naddr + 1
                      kpao = kpao + 1
                      Ssub_i(kpao,jpao) = Ssub(naddr)
                   enddo ! n3
                enddo ! n2
             endif ! End of if (label_kj_Ssub.ne.0)
          endif ! End of if (j_in_halo.le.k_in_halo)

          jpao0 = jpao0 + nd2
       enddo ! j

       kpao0 = kpao0 + nd3
    enddo ! k

    if (kpao0.ne.jpao0) then
       call cq_abort("in sub:LFD_make_Subspace_i, subspace matrix for this target atom is not square:",i)
    endif
    if (kpao0.ne.len_Sub_i) then
       call cq_abort("in sub:LFD_make_Subspace_i, size of subspace matrix for this target atom is not correct:",i)
    endif
!             
!   --- symmetrize subspace matrices ---
!
    call LFD_symm_sub_mat(Hsub_i,len_Sub_i)
    call LFD_symm_sub_mat(Ssub_i,len_Sub_i)  
!
    return
  end subroutine LFD_make_Subspace_i
!!***

!!****f* multisiteSF_module/LFD_symm_sub_mat *
!!
!!  NAME
!!   LFD_symm_sub_mat
!!
!!  PURPOSE
!!   This subroutine is for symmetrize subspace matrices,
!!   in which each value is contained only either upper or lower sides.
!!   Therefore in this subroutine, the values are copied from one side to the other side.
!!
!!   This subroutine is called in sub:LFD_make_Subspace_i
!!
!!  AUTHOR
!!   A.Nakata
!!  CREATION DATE
!!   2017/01/08
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine LFD_symm_sub_mat(mat,len)

    implicit none
    integer     :: len
    real(double) :: mat(len,len)
    integer :: i, j
    real(double) :: val

    do i = 1, len
       do j = 1, i
          val = mat(i,j) + mat(j,i)
          if (i.eq.j) val = val * half
          mat(i,j) = val
          mat(j,i) = val
       enddo
    enddo
!
    return
  end subroutine LFD_symm_sub_mat
!!***

!!****f* multisiteSF_module/LFD_make_TVEC *
!!
!!  NAME 
!!   LFD_make_TVEC
!! 
!!  PURPOSE
!!   This subroutine makes the trial vectors for this ATOM.
!!
!!   In default, trial vectors are usualy set to some PAOs on the target atom.
!!
!!   When reading from input file, the input data stored in LFD_TVEC(species,NTVEC,npao)
!!   are put into the corresponding positions of TVEC(len_Sub_i,NTVEC). 
!!
!!   This subroutine is called in sub:LocFilterDiag.
!!
!!  AUTHOR
!!   A.Nakata
!!  CREATION DATE
!!   2017/01/08
!!  MODIFICATION HISTORY
!!   2019/06/06 17:30 jack poulton and nakata
!!    introduced SZP-MSSFs with automatic setting of trial vectors using offset
!!
!!  SOURCE
!!
  subroutine LFD_make_TVEC(TVEC,NTVEC,len_Sub_i,np,i,atom_i,atom_num,atom_spec,LFDhalo, NEsub, &
                           len_Sub_i_d,n_naba_i_d,l_k_g,l_k_r2,l_kpao)

    use group_module,        ONLY: parts
    use primary_module,      ONLY: bundle
    use cover_module,        ONLY: BCS_parts
    use global_module,       ONLY: id_glob, species_glob
    use species_module,      ONLY: charge
    use matrix_data,         ONLY: mat, LD_range
    use matrix_module,       ONLY: matrix_halo
    use species_module,      ONLY: npao_species
    use pao_format
 
    implicit none

    ! Passed variables
    type(matrix_halo)  :: LFDhalo
    integer, intent(in) :: NTVEC,len_Sub_i,np,i,atom_i,atom_num,atom_spec,len_Sub_i_d,n_naba_i_d
    real(double), dimension(len_Sub_i,NTVEC) :: TVEC
    real(double) :: NEsub
    integer, dimension(n_naba_i_d) :: l_k_g
    integer, dimension(len_Sub_i_d) :: l_kpao
    real(double), dimension(n_naba_i_d) :: l_k_r2

    ! Local variables
    integer :: npao_i, &
               k, ist, gcspart, k_in_halo, npao_k, &
               neigh_global_part, atom_k, neigh_species, prncpl, prncpl0, &
               kpao, kpao0, ITVEC, ITVEC0, count, l1, nacz1, m1, k1, kk
    integer, dimension(7) :: tvec_n
    real(double) :: dx,dy,dz,r2,cutoff
    real(double), dimension(7) :: occ_n, max_cutoff_n
    real(double), parameter :: R2_ERR = 1.0e-4_double
!    real(double) :: MSSF_nonminimal_offset


    if (iprint_basis>=5.and.inode==ionode) write(io_lun,*) 'We are in sub:LFD_make_TVEC'

    npao_i = npao_species(atom_spec)

    ! --- make trial vectors
    kpao0 = 0
    kk    = 0
    ! Loop over atom_k, neighbour of atom_i
    do k = 1,mat(np,LD_range)%n_nab(i)
       ist = mat(np,LD_range)%i_acc(i) + k - 1                         ! index of atom_k in neighbour labelling
       gcspart = BCS_parts%icover_ibeg(mat(np,LD_range)%i_part(ist)) &
               + mat(np,LD_range)%i_seq(ist) -1                        ! index of atom_k in CS
       k_in_halo = LFDhalo%i_halo(gcspart)                             ! index of atom_k in halo labelling
       npao_k = LFDhalo%ndimj(k_in_halo)                               ! npao of atom_k
       
       ! Displacement vector
       dx = BCS_parts%xcover(gcspart)-bundle%xprim(atom_num)
       dy = BCS_parts%ycover(gcspart)-bundle%yprim(atom_num)
       dz = BCS_parts%zcover(gcspart)-bundle%zprim(atom_num)
       r2 = dx*dx + dy*dy + dz*dz
       r2 = sqrt(r2)
       ! check the species of neighbour
       neigh_global_part = BCS_parts%lab_cell(mat(np,LD_range)%i_part(ist))
       atom_k = id_glob(parts%icell_beg(neigh_global_part)+mat(np,LD_range)%i_seq(ist)-1)   ! global number of k
       neigh_species = species_glob(atom_k)
       ! count the total number of electrons of neighbour atoms
       NEsub = NEsub + charge(neigh_species)
       ! find on-site
       if (r2.le.R2_ERR) then
          if (atom_i.ne.atom_k) call cq_abort("in sub: make_TVEC_LFD, error in finding onsite atom", atom_i)
          if (.not.flag_LFD_ReadTVEC) then
          ! --- make trial vectors by choosing the largest PAO for each L and N
             kpao    = 0
             ITVEC   = 0
             ITVEC0  = 0
             prncpl0 = 0
             do l1 = 0, pao(atom_spec)%greatest_angmom ! Loop L
                ! check occupancy for each N to distinguish (semicore,valence) or (polarization)
                occ_n(:) = zero
                do nacz1 = 1, pao(atom_spec)%angmom(l1)%n_zeta_in_angmom
                   prncpl = pao(atom_spec)%angmom(l1)%prncpl(nacz1)
                   occ_n(prncpl) = occ_n(prncpl) + pao(atom_spec)%angmom(l1)%occ(nacz1)
                enddo
                ! find Z of the largest PAO for each L and N
                max_cutoff_n(:) = zero
                tvec_n(:) = 0
                do nacz1 = 1, pao(atom_spec)%angmom(l1)%n_zeta_in_angmom
                   prncpl = pao(atom_spec)%angmom(l1)%prncpl(nacz1)
                   if (MSSF_nonminimal_species(atom_spec).eq.1) then ! SZ-size MSSFs
                      ! only for occupied PAOs
                      if (occ_n(prncpl).gt.zero) then ! Semicore or Valence PAOs
                         cutoff = pao(atom_spec)%angmom(l1)%zeta(nacz1)%cutoff
                         if (cutoff.gt.max_cutoff_n(prncpl)) then
                            max_cutoff_n(prncpl) = cutoff
                            tvec_n(prncpl) = nacz1 ! Z of the largest PAO for this L and this N
                         endif
                      endif
                   else
                      ! for all PAOs
                      cutoff = pao(atom_spec)%angmom(l1)%zeta(nacz1)%cutoff
                      if (cutoff.gt.max_cutoff_n(prncpl)) then
                         max_cutoff_n(prncpl) = cutoff
                         tvec_n(prncpl) = nacz1 ! Z of the largest PAO for this L and this N
                      endif
                   endif
                enddo
                ! make TVEC
                do nacz1 = 1, pao(atom_spec)%angmom(l1)%n_zeta_in_angmom
                   prncpl = pao(atom_spec)%angmom(l1)%prncpl(nacz1)
                   if (nacz1.ne.1 .and. prncpl.ne.prncpl0) ITVEC0 = ITVEC0 + 2*l1 + 1
                   ! put 1 for the largest PAOs
                   if (nacz1.eq.tvec_n(prncpl)) then
                      ITVEC = ITVEC0
                      do m1 = -l1,l1
                         ITVEC = ITVEC + 1
                         kpao = kpao + 1
                         TVEC(kpao0+kpao,ITVEC) = one
                      enddo ! m1
                   else if ( MSSF_nonminimal_species(atom_spec).ne.1 .and. MSSF_nonminimal_offset.ne.zero) then
                      ITVEC = ITVEC0
                      do m1 = -l1,l1
                         ITVEC = ITVEC + 1
                         kpao = kpao + 1
                         TVEC(kpao0+kpao,ITVEC) = MSSF_nonminimal_offset
                      enddo ! m1
                   else
                      kpao = kpao + 2*l1 + 1
                   endif
                   prncpl0 = prncpl
                enddo
                ITVEC0 = ITVEC
             enddo ! l1
             if (kpao.ne.npao_k) write(io_lun,*) "in sub: make_TVEC_LFD, error in making TVEC 2", kpao, npao_k
             if (ITVEC.ne.NTVEC) write(io_lun,*) "in sub: make_TVEC_LFD, error in making TVEC 3", ITVEC, NTVEC
             if (kpao.ne.npao_k .or. ITVEC.ne.NTVEC) then
                write(io_lun,'(A/A/A)') &
                   " Most of cases, this problem occurs when the number of support function is not set correctly.",&
                   " Automatic setting of trial vectors is available only when multi-site support functions is in ",&
                   " single-zeta (SZ) size, or single-zeta plus polarization (SZP) size with flag_MSSF_nonminimal."
                call cq_abort("Stopped in sub: make_TVEC_LFD")
             endif
          else
          ! --- read TVEC from input
             do l1 =1, NTVEC
                do m1 = 1, npao_k
                   TVEC(kpao0+m1,l1) = LFD_TVEC_read(atom_spec,l1,m1)
                enddo
             enddo
          endif
       endif
       kpao0 = kpao0 + npao_k
       ! save information of neighbour atom_k to print out MOs in the subspace for debug 
       l_k_g(k)  = atom_k   ! global number of k
       l_k_r2(k) = r2       ! distance between i and k  
       do k1 = 1,npao_k
          kk = kk + 1
          l_kpao(kk) = k1
       enddo
    enddo ! k

    if (kpao0.ne.len_Sub_i) then
       call cq_abort("in sub:LFD_make_TVEC, size of subspace matrix for this target atom is not correct:",atom_i)
    endif
!
    return
  end subroutine LFD_make_TVEC
!!***

!!****f* multisiteSF_module/LFD_filter *
!!
!!  NAME 
!!   LFD_filter
!! 
!!  PURPOSE
!!   This subroutine performs filtration as
!!        t' = f(e)*t~
!!   with filter function f(e)
!!        f(e) = 1 / [exp((e-mu)/kT)+1].
!!
!!   This subroutine is called in sub:LocFilterDiag.
!!
!!  AUTHOR
!!   A.Nakata
!!  CREATION DATE
!!   2017/01/08
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine LFD_filter(TVEC,EIG,CHEMP,kT1,NUMEIG,NTVEC)                            
                            
    implicit none

    ! Passed variables
    integer, intent(in) :: NUMEIG,NTVEC
    real(double), intent(in) :: CHEMP, kT1
    real(double), intent(in), dimension(NUMEIG) :: EIG
    real(double), dimension(NUMEIG,NTVEC) :: TVEC

    ! Local variables
    real(double) :: EXFRM,FLTR
    integer :: IVEC,JVEC


    if (iprint_basis>=5.and.inode==ionode) write(io_lun,*) 'We are in sub:LFD_filter'

    do IVEC = 1,NUMEIG
       EXFRM = ( EIG(IVEC) - CHEMP ) * kT1        ! (e_i - mu) / kT
       FLTR  = one / (exp(EXFRM) + one)           ! f(e_i) = 1 / ( exp[(e_i - mu) / kT] + 1 )
       do JVEC = 1,NTVEC
          TVEC(IVEC,JVEC) = FLTR * TVEC(IVEC,JVEC)
       enddo
    enddo
!
    return
  end subroutine LFD_filter
!!***  

!!****f* multisiteSF_module/transpose_2Dmat *
!!
!!  NAME
!!   transpose_2Dmat
!!
!!  PURPOSE
!!   Transpose matrix A(i,j) into B(j,i).
!!   The size of matrices A and B should be A(LDA,LDB) and B(LDB,LDA).
!!
!!   This subroutine is called in sub:LocFilterDiag.
!!
!!  AUTHOR
!!   A.Nakata
!!  CREATION DATE
!!   2017/01/08
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine transpose_2Dmat(A,B,LDA,LDB)

    implicit none

    ! Passed variables
    integer :: LDA,LDB
    real(double) :: A(LDA,LDB), B(LDB,LDA)
    ! Passed variables    
    integer :: i,j


    if (iprint_basis>=5.and.inode==ionode) write(io_lun,*) 'We are in sub:transpose_2Dmat'
    
    do i = 1,LDA
       do j = 1,LDB
          B(j,i) = A(i,j)
       enddo
    enddo
!
    return
  end subroutine transpose_2Dmat
!!***

!!****f* multisiteSF_module/LFD_put_TVEC_to_SFcoeff *
!!
!!  NAME
!!   LFD_put_TVEC_to_SFcoeff
!!  USAGE
!!
!!  PURPOSE
!!   Put TVEC(a,l_n) into SFcoeff(i_a,k_m) if (k_m = l_n).
!!
!!   SFcoeff is in (alpha,beta,j,i) form.
!!   TVEC    is in (beta*j,alpha)   form. 
!!   Therefore, TVEC should be transposed to (alpha,beta*j) before calling this subroutine.
!!
!!   The range for Rayson's LFD should be always 
!!   equal or larger than SFcoeff_range (i.e., always dr == 0).
!!         MSmat: matrix information for Multi-site range (= SFCoeff_range)
!!        LFDmat:                        LFD_range
!!
!!   This subroutine is based on sub: addscan.
!!   The only difference is how to find the position in trial vectors,
!!   which is not in the form of CQ matrix form.
!!
!!   This subroutine is called in sub:LocFilterDiag.
!!
!!  AUTHOR
!!   A.Nakata
!!  CREATION DATE
!!   2017/01/08
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine LFD_put_TVEC_to_SFcoeff(np,i,MSmat,SFCOEFF,lenSFCOEFF,LFDmat,TVEC,lenTVEC)

    use matrix_module
    use basic_types

    implicit none

    ! Passed variables
    type(matrix) :: MSmat(:),LFDmat(:)
    integer :: lenSFCOEFF, lenTVEC
    real(double) :: SFCOEFF(lenSFCOEFF),TVEC(lenTVEC)
    integer :: np,i
    ! Local variables
    integer :: nb,ist_k,ist_l,nd1,nd2,bsize,SFCOEFFposn,LFDposn
    

    if (iprint_basis>=5.and.inode==ionode) write(io_lun,*) 'We are in sub:LFD_put_TVEC_to_SFcoeff'

    ist_l = LFDmat(np)%i_acc(i)                              ! index of atom_l in neighbor labelling for LD_range
    SFCOEFFposn = MSmat(np)%nd_offset+MSmat(np)%i_nd_acc(i)  ! starting position of matrix element for (atom_i,atom_k) in SFCOEFF matrix
    LFDposn  = 1                                             ! starting position of atom_k in TVEC matrix
    nd1 = MSmat(np)%ndimi(i)                                 ! number of SFs on atom_i (=1st dimension of matSFcoeff)

    if(MSmat(np)%n_nab(i)>0) then  ! Scan over neighbours of primary atom_i for SFcoeff_range (= atom_k)
       do nb=1,MSmat(np)%n_nab(i)
          ist_k = MSmat(np)%i_acc(i)+nb-1     ! index of atom_k in neighbour labelling for SFcoeff_range
          nd2=MSmat(np)%ndimj(ist_k)
          bsize = nd1*nd2

          ! Increment atom_l list until not greater than atom_k
          do while (LFDmat(np)%i_part(ist_l) < MSmat(np)%i_part(ist_k))
             LFDposn = LFDposn + nd1*LFDmat(np)%ndimj(ist_l)
             ist_l = ist_l+1
          enddo
          
          ! Check against the partition of atom_k
          if (LFDmat(np)%i_part(ist_l) > MSmat(np)%i_part(ist_k)) then
             cycle ! We've gone past present neighbour - increment atom_k
          else if (LFDmat(np)%i_part(ist_l)==MSmat(np)%i_part(ist_k)) then   ! atom_k and atom_l are in the same partition
             ! Check through seq and part nos
             do while (LFDmat(np)%i_seq(ist_l) < MSmat(np)%i_seq(ist_k) .AND. &
                       LFDmat(np)%i_part(ist_l) == MSmat(np)%i_part(ist_k))
                LFDposn = LFDposn + nd1*LFDmat(np)%ndimj(ist_l)
                ist_l = ist_l+1
             enddo
             ! If either partition no or seq no is too big, increment atom_k
             if (LFDmat(np)%i_part(ist_l) > MSmat(np)%i_part(ist_k)) cycle
             if (LFDmat(np)%i_seq(ist_l) > MSmat(np)%i_seq(ist_k)) then
                cycle
             else if (LFDmat(np)%i_part(ist_l)==MSmat(np)%i_part(ist_k) .AND. &
                      LFDmat(np)%i_seq(ist_l)==MSmat(np)%i_seq(ist_k)) then
                ! We've found a match - decide which matrix to store result in
                SFCOEFF(SFCOEFFposn:SFCOEFFposn+bsize-1) = TVEC(LFDposn:LFDposn+bsize-1)
                ist_l = ist_l+1
                LFDposn = LFDposn + bsize
                if (LFDposn-1>lenTVEC) call cq_abort("Length error for B in LFD_put_TVEC_to_SFcoeff: ",LFDposn-1,lenTVEC)
             endif
          endif
          SFCOEFFposn = SFCOEFFposn + bsize
          if (SFCOEFFposn-1>lenSFCOEFF) call cq_abort("Length error for A in LFD_put_TVEC_to_SFcoeff: ", &
                                                      SFCOEFFposn-1,lenSFCOEFF)
       enddo ! nb=1,MSmat(np)%n_nab(i)
    endif ! if(MSmat(np)%n_nab(i)>0)
  end subroutine LFD_put_TVEC_to_SFcoeff
!!***

!!****f* multisiteSF_module/LFD_debug_matrix *
!!
!!  NAME
!!   LFD_debug_matrix
!!  USAGE
!!
!!  PURPOSE
!!   This subroutine is used only for debug of the LFD calculation.
!!
!!   Write out a matrix (trial vectors, subspace MO coefficients, etc)
!!   not in the CQ matrix form, 
!!   with indeces 
!!   (global number, distance from atom_i, and paos for neighbour atom).
!!
!!   ITYPE : print eigenvalues or not (1/0)
!!
!!  AUTHOR
!!   A.Nakata
!!  CREATION DATE
!!   2017/01/08
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine LFD_debug_matrix(lun,itype,np,atom_i,EVAL,VEC,len_Sub_i,NUM, &
                             n_naba_i,l_k_g,l_k_r2,l_kpao)

    implicit none

    ! Local variables
    integer :: K,J,JJ,MM,MAX,JMIN,JMAX,KAT
    ! Passed variables
    integer,intent(in) :: lun, itype
    integer,intent(in) :: np, atom_i, len_Sub_i, NUM, n_naba_i
    real(double), intent(in), dimension(NUM)   :: EVAL
    real(double), intent(in), dimension(len_Sub_i,NUM) :: VEC
    integer, dimension(n_naba_i) :: l_k_g
    integer, dimension(len_Sub_i) :: l_kpao
    real(double), dimension(n_naba_i) :: l_k_r2


!   --- write out ---
    write(unit=lun,fmt='(/A,I4,5X,A,I4)') '=== atom (global): ',atom_i
       
    MAX  = 5
    JMAX = 0
    MM   = NUM / MAX + 1
    do JJ = 1, MM
       JMIN = JMAX+1
       JMAX = JMAX+MAX
       if (JMAX .gt. NUM) JMAX = NUM
       write(unit=lun,fmt='(/4x,a4,a11,a4,5(4x,i4,3x))') ' k_g','   r(i-k)  ','kpao',(J,J=JMIN,JMAX)
       if (itype.eq.1) write(unit=lun,fmt='(23x,5f11.6)') (EVAL(J),J=JMIN,JMAX)
       write(unit=lun,fmt='(A)') ' '

       KAT = 0
       do K=1,len_Sub_i
          if (l_kpao(K).eq.1) KAT = KAT + 1
          write(unit=lun,fmt='(2i4,f11.6,i4,5f11.6)') &
                K,l_k_g(KAT),l_k_r2(KAT),l_kpao(K),(VEC(K,J),J=JMIN,JMAX) 
       enddo
       if (JMAX .eq. NUM) exit
    enddo
!
    return
  end subroutine LFD_debug_matrix
!!***   

end module multisiteSF_module
