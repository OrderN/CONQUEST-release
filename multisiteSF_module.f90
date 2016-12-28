! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
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
!!   --- for setting initial coefficients ---
!!   initial_SFcoeff_SSSF, loop_initial_SFcoeff_SSSF   ! for conventional single-site SFs
!!   --- for normalise coefficients ---
!!   normalise_SFcoeff
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
!!
!!  SOURCE
!!
module multisiteSF_module

  use datatypes
  use numbers
  use global_module,          only: io_lun
  use GenComms,               only: myid, my_barrier, cq_abort, inode, ionode
  use timer_module,           only: start_timer,stop_timer
  use timer_stdclocks_module, only: tmr_std_basis,tmr_std_matrices

  implicit none

  ! RCS tag for object file identification
  character(len=80), save, private :: RCSid = "$Id$"
!!***

contains

!!****f* multisiteSF/initial_SFcoeff_SSSF *
!!
!!  NAME 
!!   initial_SFcoeff_SSSF
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
!!   This subroutine is called in sub:
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
!!
!!  SOURCE
!!
  subroutine initial_SFcoeff_SSSF
    
    use datatypes
    use numbers
    use global_module,  ONLY: iprint_basis, id_glob, species_glob, IPRINT_TIME_THRES2, nspin_SF
    use io_module,      ONLY: dump_matrix
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
    if(iprint_basis>=5.AND.myid==0) write(io_lun,'(6x,A)') ' We are in initial_SFcoeff_SSSF'


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
                   call loop_initial_SFcoeff_SSSF(iprim,halo(SFcoeff_range)%i_halo(gcspart),atom_spec)
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

    if (nspin_SF == 1) then
       call dump_matrix("SFcoeff",    matSFcoeff(1), inode)
    else
       call dump_matrix("SFcoeff_up", matSFcoeff(1), inode)
       call dump_matrix("SFcoeff_dn", matSFcoeff(2), inode)
    end if

    call stop_timer(tmr_std_matrices)
    call stop_print_timer(tmr_l_tmp1,"initial_SFcoeff_SSSF",IPRINT_TIME_THRES2)
    call stop_timer(tmr_std_basis)
!
  return
  end subroutine initial_SFcoeff_SSSF
!!***

!!****f* multisiteSF_module/loop_initial_SFcoeff_SSSF
!!
!!  NAME
!!   loop_initial_SFcoeff_SSSF
!!
!!  PURPOSE
!!   Set initial SF coefficients (SF_atom_i, PAO_atom_j)
!!   for conventional single-site SFs, in which i = j.
!!
!!   This subroutine is used instead of 
!!   sub:get_support_pao_rep in ol_rad_table_subs.f90.
!!
!!   This subroutine is called in sub:initial_SFcoeff_SSSF.
!!  
  subroutine loop_initial_SFcoeff_SSSF(iprim,j_in_halo,atom_spec)

    use global_module,  ONLY: iprint_basis, nspin_SF
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
             do nacz2 = 1, pao(atom_spec)%angmom(l2)%n_zeta_in_angmom   ! nakata3
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
  end subroutine loop_initial_SFcoeff_SSSF
!!***

!!****f* multisiteSF/normalise_SFcoeff *
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
    use global_module,  ONLY: iprint_basis, id_glob, species_glob, IPRINT_TIME_THRES2, nspin_SF
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

end module multisiteSF_module
