! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module build_PAO_matrices
! ------------------------------------------------------------------------------
! Code area 11: basis operations
! ------------------------------------------------------------------------------

!!****h* Conquest/build_PAO_matrices *
!!  NAME
!!   build_PAO_matrices
!!  PURPOSE
!!   Builds matrices from PAO integrals
!!  USES
!!
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   08:08, 2003/06/17 
!!  MODIFICATION HISTORY
!!   17:50, 2003/09/22 dave
!!    Further refinement to loop over PAO coefficients
!!  [stuck on WAGN train between New Southgate and New Barnet an hour and a half late]
!!   10:09, 13/02/2006 drb 
!!    Removed all explicit references to data_ variables and rewrote in terms of new 
!!    matrix routines
!!   2008/02/04 08:35 dave
!!    Changed for output to file not stdout
!!   2008/06/10 ast
!!    Added timers
!   2014/09/15 18:30 lat
!!    fixed call start/stop_timer to timer_module (not timer_stdlocks_module !)
!!  SOURCE
!!
module build_PAO_matrices

  use global_module,          only: io_lun
  use timer_module,           only: start_timer,stop_timer
  use timer_stdclocks_module, only: tmr_std_basis,tmr_std_matrices

  implicit none

  ! RCS tag for object file identification
  character(len=80), save, private :: RCSid = "$Id$"
!!***

contains

!!****f* build_PAO_matrices/assemble *
!!
!!  NAME 
!!   assemble
!!  USAGE
!! 
!!  PURPOSE
!!   Assembles a matrix from PAO results
!!
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   08:08, 2003/06/17 dave
!!  MODIFICATION HISTORY
!!   14:05, 2003/10/08 drb & rc
!!    Added call to calc_mat_elem and comments above
!!   16:12, 30/10/2003 drb 
!!    Checking in functioning routine.  Note that flag takes three values:
!!     1. S matrix
!!     2. KE matrix
!!     3. SP matrix (support/non-local projectors)
!!   12:54, 31/10/2003 drb 
!!    Somewhat hacked version: I've changed it so that we can do SP as well as S.  We need to sort out 
!!    how the acz loops are done, and how the coefficients are done: these are non-trivial for SP !
!!   2004/10/29 drb
!!    MAJOR bug fix: species were being done all wrong
!!   2004/11/10 drb
!!    Changed nsf association from common to maxima (as it should be !)
!!   2008/06/10 ast
!!    Added timers
!!   2009/07/08 16:48 dave
!!    Added code for one-to-one PAO to SF assignment
!!  SOURCE
!!
  subroutine assemble_2(range,matA,flag,matAD)
    
    use datatypes
    use numbers
    use global_module, ONLY: iprint_basis, id_glob, species_glob, IPRINT_TIME_THRES2
    use species_module, ONLY: species
    use GenComms, ONLY: myid, cq_abort, gcopy, my_barrier
    use matrix_module, ONLY: matrix, matrix_halo
    use group_module, ONLY: parts
    use primary_module, ONLY: bundle
    use cover_module, ONLY: BCS_parts
    use support_spec_format, ONLY: flag_paos_atoms_in_cell, flag_one_to_one
    use mult_module, ONLY: matrix_scale
    use matrix_data, ONLY: mat, halo
    use timer_module
    
    implicit none

    ! Passed variables
    integer :: length, flag, nz2, n_pfnls, i_pjnl
    integer :: range,matA
    ! Optional variables - do we want dmatrix/dPAOcoefficient ? 
    integer, OPTIONAL :: matAD

    ! Local variables
    integer :: part, memb, neigh, ist, atom_num, atom_spec, iprim, ip
    integer :: i, j, k, supfn_r, supfn_c
    integer :: neigh_global_part, neigh_global_num, neigh_species
    integer :: gcspart, wheremat, nsf_i, nsf_j, acz_i, acz_j
    integer :: l1, m1, l2, m2, i1, i2, i3, j1, j2 , j3
    real(double) :: dx, dy, dz, mat_val, tmp1, dtmp1,delt
    type(cq_timer) :: tmr_l_tmp1   

    call start_timer(tmr_std_basis)
    call start_timer(tmr_l_tmp1,WITH_LEVEL)
    if(iprint_basis>=5.AND.myid==0) write(io_lun,fmt='(6x,i5," Zeroing S")') myid
    call matrix_scale(zero,matA)
    !write(io_lun,*) 'Zeroing dS: ',PRESENT(matAD)
    if(PRESENT(matAD)) call matrix_scale(zero,matAD)
    if(iprint_basis>=5.AND.myid==0) write(io_lun,fmt='(6x,i5," Done Zeroing")') 
    iprim = 0; ip = 0
    if(flag_paos_atoms_in_cell.OR.flag==3.OR.flag_one_to_one) then
       call start_timer(tmr_std_matrices)
       do part = 1,bundle%groups_on_node ! Loop over primary set partitions
          if(iprint_basis>=6.AND.myid==0) write(io_lun,fmt='(6x,"Processor, partition: ",2i7)') myid,part
          if(bundle%nm_nodgroup(part)>0) then ! If there are atoms in partition
             do memb = 1,bundle%nm_nodgroup(part) ! Loop over atoms
                atom_num = bundle%nm_nodbeg(part)+memb-1
                iprim=iprim+1
                ip=bundle%ig_prim(iprim)
                ! Atomic species
                atom_spec = bundle%species(atom_num)
                if(iprint_basis>=6.AND.myid==0) write(io_lun,'(6x,"Processor, atom, spec: ",3i7)') myid,memb,atom_spec
                do neigh = 1, mat(part,range)%n_nab(memb) ! Loop over neighbours of atom
                   ist = mat(part,range)%i_acc(memb)+neigh-1
                   ! Build the distances between atoms - needed for phases 
                   gcspart = BCS_parts%icover_ibeg(mat(part,range)%i_part(ist))+mat(part,range)%i_seq(ist)-1
                   ! Displacement vector
                   dx = BCS_parts%xcover(gcspart)-bundle%xprim(atom_num)
                   dy = BCS_parts%ycover(gcspart)-bundle%yprim(atom_num)
                   dz = BCS_parts%zcover(gcspart)-bundle%zprim(atom_num)
                   if(iprint_basis>=6.AND.myid==0) write(80+myid,*) 'dx,y,z: ',gcspart,dx,dy,dz
                   ! We need to know the species of neighbour
                   neigh_global_part = BCS_parts%lab_cell(mat(part,range)%i_part(ist)) 
                   neigh_global_num  = id_glob(parts%icell_beg(neigh_global_part)+mat(part,range)%i_seq(ist)-1)
                   neigh_species = species_glob(neigh_global_num)
                   if(flag_one_to_one) then
                      ip = atom_spec
                      neigh_global_num = neigh_species
                   end if
                   if(iprint_basis>=6.AND.myid==0) write(io_lun,'(6x,"Processor, neighbour, spec: ",3i7)') myid,neigh,neigh_species
                   ! Now loop over support functions and PAOs and call routine
                   if(flag<3) then
                      if(PRESENT(matAD)) then
                         call loop_12(matA,iprim,halo(range)%i_halo(gcspart),ip,neigh_global_num,&
                              dx,dy,dz,flag ,atom_spec,neigh_species,0,0,matAD)
                      else
                         call loop_12(matA,iprim,halo(range)%i_halo(gcspart),ip,neigh_global_num,&
                              dx,dy,dz,flag ,atom_spec,neigh_species,0,0)
                      end if
                   else
                      if(PRESENT(matAD)) then
                         call loop_3(matA,iprim,halo(range)%i_halo(gcspart),ip,neigh_global_num,dx,dy,dz,&
                              &atom_spec,neigh_species,0,0,matAD)
                      else
                         call loop_3(matA,iprim,halo(range)%i_halo(gcspart),ip,neigh_global_num,dx,dy,dz,&
                              &atom_spec,neigh_species,0,0)
                      end if
                   endif
                end do ! End do neigh=1,mat%n_nab
             end do ! End do memb =1,nm_nodgroup
          end if ! End if nm_nodgroup > 0
       end do ! End do part=1,groups_on_node
       call stop_timer(tmr_std_matrices)
    else
       call cq_abort("Storage of PAOs on primary-set processors not implemented yet")
    end if    
    call stop_print_timer(tmr_l_tmp1,"S matrix assembly",IPRINT_TIME_THRES2)
    call stop_timer(tmr_std_basis)
  end subroutine assemble_2
!!***
  
!!   2009/07/08 16:48 dave
!!    Added code for one-to-one PAO to SF assignment
  subroutine loop_12(matA,iprim,j_in_halo,atom_i,atom_j,dx,dy,dz,flag&
       &,atom_spec,neigh_species,deriv_flag,direction,matAD)

    use datatypes
    use support_spec_format, ONLY: supports_on_atom, flag_paos_atoms_in_cell, supports_on_atom_remote, &
         support_function, flag_one_to_one
    use angular_coeff_routines, ONLY: calc_mat_elem_gen, grad_mat_elem_gen2_prot
    use GenComms, ONLY: myid, cq_abort
    use mult_module, ONLY: store_matrix_value_pos, matrix_pos, return_matrix_value_pos

    implicit none

    !routine to calculate the case(1) pao_pao and case(2) pao_ke_pao
    !matrix elements
    integer, intent(in) ::  iprim,j_in_halo,atom_j,atom_i,flag,matA
    integer, intent(in) :: atom_spec,neigh_species,deriv_flag,direction
    real(double), intent(in) :: dx,dy,dz
    ! Optional variables - do we want dmatrix/dPAOcoefficient ? 
    integer, OPTIONAL :: matAD
    ! Local variables
    integer :: wheremat
    real(double) :: mat_val,mat_val_dummy, coeff_j
    integer :: i1,j1,k1,l1,m1,i2,j2,i3,j3,k2,l2,m2,count1,count2,nacz1,nacz2
    integer :: i4,j4, n_support1,n_support2,myflag
    type(support_function), pointer, dimension(:) :: supports_on_atom_j

    mat_val = 0.0_double
    count1 = 1
    if(flag_paos_atoms_in_cell) then
       supports_on_atom_j => supports_on_atom
    else
       supports_on_atom_j => supports_on_atom_remote
    end if
    ! Loop over PAOs on i
    do l1 = 0, supports_on_atom(atom_i)%lmax
       do nacz1 = 1, supports_on_atom(atom_i)%naczs(l1)
          do m1 = -l1,l1
             count2 = 1
             ! Loop over PAOs on j
             do l2 = 0, supports_on_atom_j(atom_j)%lmax
                do nacz2 = 1, supports_on_atom_j(atom_j)%naczs(l2)
                   do m2 = -l2,l2
                      if(deriv_flag.eq.0) then
                         call calc_mat_elem_gen(flag,atom_spec,l1,nacz1,m1,neigh_species, &
                              l2,nacz2,m2,dx,dy,dz,mat_val)
                      else if(deriv_flag.eq.1) then
                         call grad_mat_elem_gen2_prot(direction,flag,atom_spec,l1,nacz1,m1,neigh_species, &
                              l2,nacz2,m2,dx,dy,dz,mat_val)
                      endif
                      if(flag_one_to_one) then
                         mat_val_dummy = mat_val
                         wheremat = matrix_pos(matA,iprim,j_in_halo,count1,count2)
                         call store_matrix_value_pos(matA,wheremat,mat_val)
                         if(PRESENT(matAD).AND.count1==1) then
                            wheremat = matrix_pos(matAD,iprim,j_in_halo,count1,count2)
                            call store_matrix_value_pos(matAD,wheremat,mat_val)

                         end if
                      else
                         do n_support1 = 1, supports_on_atom(atom_i)%nsuppfuncs
                            do n_support2 = 1, supports_on_atom_j(atom_j)%nsuppfuncs
                               !--------------------------------------------------------!
                               !     Now multiply by corresponding pao coefficients     !
                               !--------------------------------------------------------!
                               coeff_j = supports_on_atom_j(atom_j)%supp_func(n_support2)%coefficients(count2)
                               mat_val_dummy = mat_val* &
                                    supports_on_atom(atom_i)%supp_func(n_support1)%coefficients(count1)*coeff_j
                               wheremat = matrix_pos(matA,iprim,j_in_halo,n_support1,n_support2)
                               call store_matrix_value_pos(matA,wheremat,mat_val_dummy)
                               if(PRESENT(matAD).AND.n_support1==1) then
                                  mat_val_dummy = mat_val*coeff_j
                                  wheremat = matrix_pos(matAD,iprim,j_in_halo,count1,n_support2)
                                  call store_matrix_value_pos(matAD,wheremat,mat_val_dummy)
                               end if
                            enddo ! nsupport2
                         enddo ! nsupport1
                      end if
                      count2 = count2 +1
                   enddo ! m2
                enddo ! nacz2
             enddo !l2
             count1 = count1 + 1
          enddo ! m1
       enddo ! nacz1
    enddo ! l1
    nullify(supports_on_atom_j)
  end subroutine loop_12
  
!!   2009/07/08 16:48 dave
!!    Added code for one-to-one PAO to SF assignment
  subroutine loop_3(matA,iprim,j_in_halo,atom_i,atom_j,dx,dy,dz&
       &,atom_spec,neigh_species,deriv_flag,direction,matAD)

    use datatypes
    use pseudo_tm_info, ONLY: pseudo
    use support_spec_format, ONLY: supports_on_atom, flag_one_to_one !contains all info on support-pao representation
    use angular_coeff_routines, ONLY: calc_mat_elem_gen, grad_mat_elem_gen2_prot
    use GenComms, ONLY: myid, cq_abort
    use mult_module, ONLY: store_matrix_value_pos, matrix_pos, return_matrix_value_pos

    implicit none

    !routine to generate <nlpf|pao> matrix elements
    integer, intent(in) ::  iprim,j_in_halo,atom_i,atom_j, matA
    integer, intent(in) :: atom_spec,neigh_species,deriv_flag,direction
    real(double), intent(in) :: dx,dy,dz
    ! Optional variables - do we want dmatrix/dPAOcoefficient ? 
    integer, OPTIONAL :: matAD
    ! Local variables
    real(double) :: mat_val,mat_val_dummy
    integer :: i1,k1,l1,m1,i2,i3,i4,k2,l2,m2,count1,count2,nacz1,nacz2,n,n_support1
    integer :: n_pfnls,i_pfnl,nzeta2,index_j,index_i,part_index_j,wheremat

    !now the loop over support functions using the new datatype

    n_pfnls = pseudo(neigh_species)%n_pjnl !no of projector functions
    mat_val = 0.0_double

    count1 = 1
    do l1 = 0,supports_on_atom(atom_i)%lmax
       do nacz1 = 1, supports_on_atom(atom_i)%naczs(l1)
          do m1 = -l1,l1
             ! Fix DRB 2007/03/22
             part_index_j = 0
            if(n_pfnls>0) then
             do i_pfnl = 1, n_pfnls
                l2 = pseudo(neigh_species)%pjnl_l(i_pfnl)
                nzeta2 = pseudo(neigh_species)%pjnl_n(i_pfnl)
                do m2 = -l2,l2
                   if(deriv_flag.eq.0) then
                      call calc_mat_elem_gen(3,atom_spec,l1,nacz1,m1,neigh_species,l2,nzeta2,m2,dx,dy,dz,mat_val)
                   else if(deriv_flag.eq.1) then
                      call grad_mat_elem_gen2_prot(direction,3,atom_spec,l1,nacz1,m1,neigh_species,&
                           l2,nzeta2,m2,dx,dy,dz,mat_val)
                   endif
                   ! Fix DRB 2007/03/22
                   !part_index_j = 0
                   !if(l2>0) then
                   !   do n = 0,l2-1
                   !      part_index_j = part_index_j+(2*n+1)
                   !   enddo
                   !end if
                   index_j = part_index_j + (m2+l2+1)
                   if(supports_on_atom(atom_i)%nsuppfuncs > 0) then
                      if(flag_one_to_one) then
                         wheremat = matrix_pos(matA,iprim,j_in_halo,count1,index_j)
                         call store_matrix_value_pos(matA,wheremat,mat_val)
                      else
                         do n_support1 = 1, supports_on_atom(atom_i)%nsuppfuncs
                            mat_val_dummy = mat_val*supports_on_atom(atom_i)%supp_func(n_support1)%coefficients(count1)
                            wheremat = matrix_pos(matA,iprim,j_in_halo,n_support1,index_j)
                            call store_matrix_value_pos(matA,wheremat,mat_val_dummy)
                         enddo ! n_support1
                      end if
                   endif
                   if(PRESENT(matAD)) then
                      wheremat = matrix_pos(matAD,iprim,j_in_halo,count1,index_j)
                      call store_matrix_value_pos(matAD,wheremat,mat_val)
                   end if
                enddo ! m2
                ! Fix DRB 2007/03/22
                part_index_j = part_index_j + (2*l2+1)
             enddo ! i_pfnl
             count1 = count1+1
            endif ! (n_pfnls>0) then
          enddo ! m1
       enddo ! nacz1
    enddo ! l1
  end subroutine loop_3
    


!!****f* build_PAO_matrices/assemble_deriv *
!!
!!  NAME 
!!   assemble_deriv
!!  USAGE
!! 
!!  PURPOSE
!!   Assembles derivative of matrix from PAO results
!!
!!  INPUTS
!! 
!! 
!!  USES
!! 
!!  AUTHOR
!!   D.R.Bowler
!!  CREATION DATE
!!   17:17, 2003/12/02 
!!  MODIFICATION HISTORY
!!   2004/10/29 drb
!!    MAJOR bug fix: species were being done all wrong, and wrong gradient being called
!!   2004/11/10 drb
!!    Changed nsf association from common to maxima (as it should be !)
!!   2008/06/10 ast
!!    Added timers
!!   2009/07/08 16:48 dave
!!    Added code for one-to-one PAO to SF assignment
!!  SOURCE
!!
  subroutine assemble_deriv_2(direction,range,matA,flag)
    
    use datatypes
    use numbers
    use global_module, ONLY: iprint_basis, id_glob, species_glob
    use species_module, ONLY: species
    use GenComms, ONLY: myid, cq_abort, gcopy
    use matrix_module, ONLY: matrix, matrix_halo
    use group_module, ONLY: parts
    use primary_module, ONLY: bundle
    use cover_module, ONLY: BCS_parts
    use angular_coeff_routines, ONLY: grad_mat_elem_gen2_prot
    !use angular_coeff_routines, ONLY: numerical_ol_gradient
    use matrix_data, ONLY: mat, halo
    use mult_module, ONLY: matrix_scale, store_matrix_value_pos, matrix_pos
    use support_spec_format, ONLY: flag_paos_atoms_in_cell, flag_one_to_one
    
    implicit none

    ! Passed variables
    integer :: direction, flag, range, matA

    ! Local variables
    integer :: part, memb, neigh, ist, atom_num, atom_spec, iprim, ip
    integer :: i, j, k, supfn_r, supfn_c
    integer :: neigh_global_part, neigh_global_num, neigh_species
    integer :: gcspart, wheremat, nsf_i, nsf_j, acz_i, acz_j
    integer :: l1, m1, l2, m2, n1, n2
    real(double) :: dx, dy, dz, mat_val, mat_val2, mat_val3, dr
    
    logical :: paoloop ! Use new loop or old loop ?
    
    call start_timer(tmr_std_matrices)
    call start_timer(tmr_std_basis)
    call matrix_scale(zero,matA)
    iprim = 0; ip = 0
    if(flag_paos_atoms_in_cell.OR.flag==3.OR.flag_one_to_one) then
       do part = 1,bundle%groups_on_node ! Loop over primary set partitions
          if(iprint_basis>=6.AND.myid==0) write(io_lun,fmt='(6x,"Processor, partition: ",2i7)') myid,part
          if(bundle%nm_nodgroup(part)>0) then ! If there are atoms in partition
             do memb = 1,bundle%nm_nodgroup(part) ! Loop over atoms
                atom_num = bundle%nm_nodbeg(part)+memb-1
                !RC setting up so i have the global number of the bra atom  
                iprim=iprim+1
                ip=bundle%ig_prim(iprim)
                ! Atomic species
                atom_spec = bundle%species(atom_num)
                if(iprint_basis>=6.AND.myid==0) write(io_lun,'(6x,"Processor, atom, spec: ",3i7)') myid,memb,atom_spec
                do neigh = 1, mat(part,range)%n_nab(memb) ! Loop over neighbours of atom
                   ist = mat(part,range)%i_acc(memb)+neigh-1
                   ! Build the distances between atoms - needed for phases 
                   gcspart = BCS_parts%icover_ibeg(mat(part,range)%i_part(ist))+mat(part,range)%i_seq(ist)-1
                   ! Displacement vector
                   dx = BCS_parts%xcover(gcspart)-bundle%xprim(atom_num)
                   dy = BCS_parts%ycover(gcspart)-bundle%yprim(atom_num)
                   dz = BCS_parts%zcover(gcspart)-bundle%zprim(atom_num)
                   if(iprint_basis>=6.AND.myid==0) write(80+myid,*) 'dx,y,z: ',gcspart,dx,dy,dz
                   ! We need to know the species of neighbour
                   neigh_global_part = BCS_parts%lab_cell(mat(part,range)%i_part(ist)) 
                   neigh_global_num  = id_glob(parts%icell_beg(neigh_global_part)+mat(part,range)%i_seq(ist)-1)
                   neigh_species = species_glob(neigh_global_num)
                   if(flag_one_to_one) then
                      ip = atom_spec
                      neigh_global_num = neigh_species
                   end if
                   ! Where to put the result - this will become matrix_pos
                   wheremat = matrix_pos(matA,iprim,halo(range)%i_halo(gcspart))
                   ! Now loop over support functions and PAOs and call routine
                   if(wheremat/=mat(part,range)%onsite(memb)) then
                      if(flag<3) then
                         call loop_12(matA,iprim,halo(range)%i_halo(gcspart),ip,neigh_global_num,dx,dy,dz,flag&
                              &,atom_spec,neigh_species,1,direction)
                      else
                         call loop_3(matA,iprim,halo(range)%i_halo(gcspart),ip,neigh_global_num,dx,dy,dz,&
                              &atom_spec,neigh_species,1,direction)
                      endif
                   end if
                end do ! End do neigh=1,mat%n_nab
             end do ! End do memb =1,nm_nodgroup
          end if ! End if nm_nodgroup > 0
       end do ! End do part=1,groups_on_node
    else
       call cq_abort("Storage of PAOs on primary-set processors not implemented yet")
    end if
    call stop_timer(tmr_std_basis)
    call stop_timer(tmr_std_matrices)
  end subroutine assemble_deriv_2
!!***

!%%! ! Temp for building s and t for preconditioning
!%%! 
!%%!   subroutine assemble_pao(mat,data_build,length,nf1,nf2,flag)
!%%!     
!%%!     use datatypes
!%%!     use global_module, ONLY: iprint_basis, id_glob, species_glob
!%%!     use species_module, ONLY: species
!%%!     use GenComms, ONLY: myid, cq_abort, gcopy, my_barrier
!%%!     use pseudo_tm_info, ONLY: pseudo
!%%!     use common, ONLY: nsf
!%%!     use maxima_module, ONLY: mx_nponn
!%%!     use matrix_module, ONLY: matrix
!%%!     use group_module, ONLY: parts
!%%!     use primary_module, ONLY: bundle
!%%!     use cover_module, ONLY: BCS_parts
!%%!     use angular_coeff_routines, ONLY: calc_mat_elem_gen
!%%!     !use support_spec_format, ONLY : which_support
!%%!     use support_spec_format, ONLY: supports_on_atom
!%%!     
!%%!     implicit none
!%%! 
!%%!     ! Passed variables
!%%!     integer :: length, flag, nf1, nf2, nz2, n_pfnls, i_pjnl
!%%!     type(matrix), dimension(mx_nponn), intent(in) :: mat
!%%!     real(double) :: data_build(nf1,nf2,length)
!%%! 
!%%!     ! Local variables
!%%!     integer :: part, memb, neigh, ist, atom_num, atom_spec, iprim, ip
!%%!     integer :: i, j, k, supfn_r, supfn_c
!%%!     integer :: neigh_global_part, neigh_global_num, neigh_species
!%%!     integer :: gcspart, wheremat, nsf_i, nsf_j, acz_i, acz_j
!%%!     integer :: l1, m1, l2, m2, i1, i2, i3, j1, j2 , j3
!%%!     real(double) :: dx, dy, dz, mat_val, tmp1, dtmp1,delt
!%%!     
!%%!     data_build = 0.0_double
!%%!     iprim = 0; ip = 0
!%%!     do part = 1,bundle%groups_on_node ! Loop over primary set partitions
!%%!        if(iprint_basis>=6) write(io_lun,fmt='(6x,"Processor, partition: ",2i7)') myid,part
!%%!        if(bundle%nm_nodgroup(part)>0) then ! If there are atoms in partition
!%%!           do memb = 1,bundle%nm_nodgroup(part) ! Loop over atoms
!%%!              atom_num = bundle%nm_nodbeg(part)+memb-1
!%%!              iprim=iprim+1
!%%!              ip=bundle%ig_prim(iprim)
!%%!              ! Atomic species
!%%!              atom_spec = bundle%species(atom_num)
!%%!              if(iprint_basis>=6) write(io_lun,'(6x,"Processor, atom, spec: ",3i7)') myid,memb,atom_spec
!%%!              do neigh = 1, mat(part)%n_nab(memb) ! Loop over neighbours of atom
!%%!                 ist = mat(part)%i_acc(memb)+neigh-1
!%%!                 ! Build the distances between atoms - needed for phases 
!%%!                 gcspart = BCS_parts%icover_ibeg(mat(part)%i_part(ist))+mat(part)%i_seq(ist)-1
!%%!                 ! Displacement vector
!%%!                 dx = BCS_parts%xcover(gcspart)-bundle%xprim(atom_num)
!%%!                 dy = BCS_parts%ycover(gcspart)-bundle%yprim(atom_num)
!%%!                 dz = BCS_parts%zcover(gcspart)-bundle%zprim(atom_num)
!%%!                 if(iprint_basis>=6) write(80+myid,*) 'dx,y,z: ',gcspart,dx,dy,dz
!%%!                 ! We need to know the species of neighbour
!%%!                 neigh_global_part = BCS_parts%lab_cell(mat(part)%i_part(ist)) 
!%%!                 neigh_global_num  = id_glob(parts%icell_beg(neigh_global_part)+mat(part)%i_seq(ist)-1)
!%%!                 neigh_species = species_glob(neigh_global_num)
!%%!                 if(iprint_basis>=6) write(io_lun,'(6x,"Processor, neighbour, spec: ",3i7)') myid,neigh,neigh_species
!%%!                 ! Where to put the result
!%%!                 wheremat = mat(part)%offset + ist
!%%!                 if(wheremat>mat(part)%length) call cq_abort("Error in matrix length ",mat(part)%length,wheremat)
!%%!                 ! Now loop over support functions and PAOs and call routine
!%%!                 call loop_pao(data_build,nf1,nf2,length,wheremat,ip,neigh_global_num,dx,dy,dz,flag ,atom_spec,neigh_species)
!%%!              end do ! End do neigh=1,mat%n_nab
!%%!           end do ! End do memb =1,nm_nodgroup
!%%!        end if ! End if nm_nodgroup > 0
!%%!     end do ! End do part=1,groups_on_node
!%%!     !RC putting in a stop sign
!%%!     
!%%!   end subroutine assemble_pao
!%%!   
!%%!   subroutine loop_pao(data_build,nf1,nf2,length,wheremat,global_i,global_j,dx,dy,dz,flag ,atom_spec,neigh_species)
!%%!     use datatypes
!%%!     use support_spec_format, ONLY: which_support, supports_on_atom
!%%!     use angular_coeff_routines, ONLY: calc_mat_elem_gen, grad_mat_elem_gen2_prot
!%%!     use GenComms, ONLY: myid, cq_abort
!%%!     implicit none
!%%!     !routine to calculate the case(1) pao_pao and case(2) pao_ke_pao
!%%!     !matrix elements
!%%!     integer, intent(in) ::  wheremat,global_j,global_i,nf1,nf2,length,flag
!%%!     integer, intent(in) :: atom_spec,neigh_species
!%%!     real(double), intent(in) :: dx,dy,dz
!%%!     real(double) :: data_build(nf1,nf2,length)
!%%!     real(double) :: mat_val,mat_val_dummy
!%%!     integer :: i1,j1,k1,l1,m1,i2,j2,i3,j3,k2,l2,m2,count1,count2,nacz1,nacz2
!%%!     integer :: i4,j4, n_support1,n_support2,myflag
!%%!     
!%%!     mat_val = 0.0_double
!%%!     count1 = 1
!%%!     ! Loop over PAOs on i
!%%!     do l1 = 0, supports_on_atom(global_i)%lmax
!%%!        do nacz1 = 1, supports_on_atom(global_i)%naczs(l1)
!%%!           do m1 = -l1,l1
!%%!              count2 = 1
!%%!              ! Loop over PAOs on j
!%%!              do l2 = 0, supports_on_atom(global_j)%lmax
!%%!                 do nacz2 = 1, supports_on_atom(global_j)%naczs(l2)
!%%!                    do m2 = -l2,l2
!%%!                       call calc_mat_elem_gen(flag,atom_spec,l1,nacz1,m1,neigh_species, &
!%%!                            l2,nacz2,m2,dx,dy,dz,mat_val)
!%%!                       data_build(count1,count2,wheremat) = data_build(count1,count2,wheremat) + mat_val
!%%!                       count2 = count2 +1
!%%!                    enddo ! m2
!%%!                 enddo ! nacz2
!%%!              enddo !l2
!%%!              count1 = count1 + 1
!%%!           enddo ! m1
!%%!        enddo ! nacz1
!%%!     enddo ! l1
!%%!   end subroutine loop_pao
!%%!   
  

end module build_PAO_matrices
