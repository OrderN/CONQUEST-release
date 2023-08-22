! -*- mode: F90; mode: font-lock -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module mlff
! ------------------------------------------------------------------------------
! (Code area 9: general) not sure ***
! ------------------------------------------------------------------------------

!!****h* Conquest/mlff
!!  NAME
!!   mlff
!!  PURPOSE
!!   Obtain forces and energies from machine learning force field
!!  AUTHOR
!!   Jianbo.Lin
!!  CREATION DATE
!!   2022/03/14
!!  MODIFICATION HISTORY
!
!!  SOURCE
!!
module mlff

  use datatypes
  use energy, only: ml_energy

  implicit none

  save

  real(double), allocatable :: ml_force(:, :)
  real(double), dimension(3,3) :: ml_stress


contains


!!****f* mlff_module/matrix_ini_ML
!! modified from matrix_elements_module/matrix_ini
!!  NAME
!!   matrix_ini_ML
!!  USAGE
!!
!!  PURPOSE
!!   Creates indexing for matrix
!!  INPUTS
!!
!!
!!  USES
!!
!!  AUTHOR
!!   JianBo.Lin
!!  CREATION DATE
!!   25/11/99
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
!  subroutine matrix_ini_ML(parts,prim,gcs,amat,amat_ML,rcut,myid)
  subroutine matrix_ini_ML(parts,prim,gcs,amat,amat_ML, rcut,myid)

    ! Module usage
    use datatypes
    use global_module
    use basic_types
    use matrix_module
    use trans_module
    use GenComms, ONLY: my_barrier, gmax
    use mlff_type, ONLY: matrix_ML, allocate_matrix_ML,params_ML,&
            allocate_matrix_ML_triplet
    use matrix_elements_module, ONLY: make_npxyz, get_naba_max

    implicit none

    ! Passed variables
    ! First, the small group variables
    type(group_set) :: parts
    type(primary_set) :: prim
    type(cover_set) :: gcs
    ! Now general integers
    integer :: myid
    real(double) :: rcut
    ! These variables are generic matrix variables
    type(matrix), dimension (:) :: amat
    type(matrix_ML), dimension (:) :: amat_ML
    !!integer(integ), pointer, dimension(:) :: aind

    ! Local variables
    integer :: nnd
    integer :: ierr
    integer :: index_size
    integer, dimension(:), allocatable :: part_offset

    nnd = myid+1
    if(prim%n_prim.gt.0) then ! If we actually HAVE a primary set
       if(.not.allocated(part_offset)) allocate(part_offset(prim%groups_on_node),STAT=ierr)
       if(ierr/=0) call cq_abort("Error allocating part_offset in matrix_ini_ML: ", prim%groups_on_node,ierr)
       call get_naba_max(prim,gcs,amat,rcut,index_size,part_offset)
       ! The matrix maxima MUST be global across processors !
       call gmax(amat(1)%mx_nab)
       amat(2:prim%groups_on_node)%mx_nab = amat(1)%mx_nab
       amat_ML%mx_nab = amat(1)%mx_nab
       call gmax(amat(1)%mx_abs)
       amat(2:prim%groups_on_node)%mx_abs = amat(1)%mx_abs
       amat_ML%mx_abs = amat(1)%mx_abs
       !write(io_lun,*) 'Matrix maxima: ',amat(1)%mx_nab, amat(1)%mx_abs
       !write(io_lun,*) 'Index length: ',parts%mx_ngonn*(3*parts%mx_mem_grp+5*parts%mx_mem_grp*amat(1)%mx_abs)
       !allocate(aind(parts%mx_ngonn*(3*parts%mx_mem_grp+5*parts%mx_mem_grp*amat(1)%mx_abs)),STAT=ierr)
       !if(iprint_index > 1) write(io_lun,*) 'Matrix index: Processor ',myid,'Size ',index_size
       !!allocate(aind(index_size),STAT=ierr)
       !!if(ierr/=0) call cq_abort("Error allocating matrix index in matrix_ini_ML: ", index_size,ierr)
            !parts%mx_ngonn*(3*parts%mx_mem_grp+5*parts%mx_mem_grp*amat(1)%mx_abs),ierr)
       !!aind = 0
       !!call set_matrix_pointers(prim%nm_nodgroup,amat, aind, prim%groups_on_node,parts%mx_mem_grp,part_offset)
       !!call allocate_matrix(amat, prim%groups_on_node,parts%mx_mem_grp)
       call allocate_matrix_ML(amat_ML, prim%groups_on_node,parts%mx_mem_grp)
       call get_naba_ML(prim,gcs,amat_ML,rcut)

       if(index(params_ML%descriptor_type, '3b') .gt. 0)then
          write(*,*) '3b is found at descriptor_type: ', &
                  index(params_ML%descriptor_type, '3b')
          call allocate_matrix_ML_triplet(amat_ML, prim%groups_on_node)
          call get_triplet_ML(amat_ML, prim%groups_on_node, rcut)
       end if

       ! Now that we know neighbour numbers, sort out the pointers
       !!call set_matrix_pointers2(prim%nm_nodgroup,amat,aind, &
       !!     prim%groups_on_node,parts%mx_mem_grp)
       ! * IMPORTANT * ! This must come next: it builds ndimj
       !!call make_npxyz(nnd,amat,prim,gcs,parts%mx_mem_grp)
       call my_barrier
    endif ! n_prim.gt.0
    return
  end subroutine matrix_ini_ML
!!***



!!****f*  ML_type/read_split *
!!
!!  NAME
!!   read_split
!!  PURPOSE
!!
!!  AUTHOR
!!   Jianbo Lin
!!  CREATION DATE
!!   2022/03/14
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!


!!****f* mlff_module/get_naba_ML *
!! modified from matrix_elements_module/get_naba
!!
!!  NAME
!!   get_naba_ML
!!  USAGE
!!
!!  PURPOSE
!!   Creates neighbour lists for a matrix range rcut
!!  INPUTS
!!   type(primary_set) :: prim
!!   type(cover_set) :: gcs
!!   type(matrix), dimension(:) :: amat
!!   real(double) :: rcut
!!  USES
!!
!!  AUTHOR
!!   JianBo.Lin
!!  CREATION DATE
!!   2022/07/20
!!  MODIFICATION HISTORY
!!   2022/07/27 JL
!!    Modified to get pair information for machine learing
!!  SOURCE
!!
  subroutine get_naba_ML(prim,gcs,amat,rcut)

    ! Module usage
    use datatypes
    use basic_types
    use matrix_module
    use GenComms, ONLY: cq_abort, inode, ionode, my_barrier
    use global_module, ONLY: id_glob, species_glob, sf, nlpf, paof, napf, numprocs
    use group_module, ONLY: parts
    use species_module, ONLY: nsf_species, nlpf_species, npao_species, napf_species
    use mlff_type

    implicit none

    ! Passed variables
    type(primary_set), intent(in) :: prim
    type(cover_set), intent(in) :: gcs
    real(double), intent(in) :: rcut
    type(matrix_ML), dimension (:), intent(inout) :: amat

    ! Local variables
    integer :: inp, cumu_ndims, neigh_spec, ip_process
    integer :: nn,j,np,ni,ist
    real(double) :: rcutsq,dx,dy,dz,rij2
    real(double), parameter :: tol=1.0e-8_double

    ! Check that prim and gcs are correctly set up
    if((.NOT.ASSOCIATED(gcs%xcover)).OR. &
         (.NOT.ASSOCIATED(prim%xprim))) then
       call cq_abort('get_naba: gcs or prim without members !')
    endif
    rcutsq=rcut*rcut

    ! loop over all atom pairs (atoms in primary set, max. cover set) -
    inp=1  ! Indexes primary atoms
    amat(1)%offset=0
    amat(1)%nd_offset=0
    do nn=1,prim%groups_on_node ! Partitions in primary set
       !check get_naba_ML!write(*,*) 'nn loop, inode=',inode

       amat(nn)%part_nabs = 0    ! Cumulative neighbours of partition
       amat(nn)%part_nd_nabs = 0    ! Cumulative neighbours of partition
       if(prim%nm_nodgroup(nn).gt.0) then  ! Are there atoms ?
          amat(nn)%i_acc(1)=1
          amat(nn)%i_nd_acc(1)=1
          amat(nn)%n_atoms = prim%nm_nodgroup(nn) ! Redundant, but useful
          !write(io_lun,*) 'Starting group with atoms: ',nn,prim%nm_nodgroup(nn)
          do j=1,prim%nm_nodgroup(nn)  ! Loop over atoms in partition
             amat(nn)%n_nab(j)=0
             amat(nn)%i_species(j)=prim%species(inp)

             cumu_ndims = 0
             do np=1,gcs%ng_cover  ! Loop over partitions in GCS
                if(gcs%n_ing_cover(np).gt.0) then  ! Are there atoms ?
                   !if(gcs%icover_ibeg(np)+gcs%n_ing_cover(np)-1.gt.&
                   !     gcs%mx_mcover) then
                   !   call cq_abort('get_naba: overran gcs mx_mcover: ', &
                   !        gcs%icover_ibeg(np)+gcs%n_ing_cover(np)-1, &
                   !        gcs%mx_mcover)
                   !endif
                   do ni=1,gcs%n_ing_cover(np)
                      dx=gcs%xcover(gcs%icover_ibeg(np)+ni-1)-prim%xprim(inp)
                      dy=gcs%ycover(gcs%icover_ibeg(np)+ni-1)-prim%yprim(inp)
                      dz=gcs%zcover(gcs%icover_ibeg(np)+ni-1)-prim%zprim(inp)
                      rij2=dx*dx+dy*dy+dz*dz
                      if(rij2<rcutsq-tol) then
                         !write(*,*) 'gcs%lab_cover(np),gcs%iprim_group(inp), inode=',inode,&
                         !        gcs%lab_cover(np), gcs%iprim_group(inp)
                         if(gcs%lab_cover(np)==gcs%iprim_group(inp).AND.ni==j) then
                            amat(nn)%onsite(j)=amat(nn)%nd_offset+amat(nn)%i_nd_acc(j)+cumu_ndims
                            cycle                                   ! no need to collect onsite atoms for ML
                         endif

                         amat(nn)%n_nab(j)=amat(nn)%n_nab(j)+1
                         !write(io_lun,*) 'Neighbour: ',j,amat(nn)%n_nab(j)
                         if(amat(nn)%n_nab(j).gt.amat(nn)%mx_abs) then
                            call cq_abort('get_naba: n_nab>mx_nab: ',amat(nn)%n_nab(j), &
                                 amat(nn)%mx_abs)
                         endif
                         ist = amat(nn)%i_acc(j)+amat(nn)%n_nab(j)-1
                         amat(nn)%i_part(ist)=np
                         amat(nn)%i_seq(ist)=ni
                         amat(nn)%radius(ist)=sqrt(rij2)
                         amat(nn)%dx(ist)=dx
                         amat(nn)%dy(ist)=dy
                         amat(nn)%dz(ist)=dz

                         neigh_spec = species_glob( id_glob( parts%icell_beg(gcs%lab_cell(np)) +ni-1 ))
                         amat(nn)%j_species(ist)=neigh_spec

                      endif
                   enddo ! End n_inp_cover
                endif
             enddo ! End np_cover
             !write(io_lun,*) 'Finishing prim atom: ',inode,inp,cumu_ndims,j,prim%nm_nodgroup(nn)
             !check get_naba_ML!write(*,*) 'after np_cover', inode
             if(j.lt.prim%nm_nodgroup(nn)) then
                amat(nn)%i_acc(j+1)=amat(nn)%i_acc(j)+amat(nn)%n_nab(j)
                amat(nn)%i_nd_acc(j+1)=amat(nn)%i_nd_acc(j)+cumu_ndims
             endif
             amat(nn)%part_nabs = amat(nn)%part_nabs+amat(nn)%n_nab(j)
             amat(nn)%part_nd_nabs = amat(nn)%part_nd_nabs+cumu_ndims
             inp=inp+1  ! Indexes primary-set atoms
             !write(*,*) 'after part_nd_nabs in get_naba_ML', inode
          enddo ! End prim%nm_nodgroup
          !write(*,*) 'after prim%nm_nodgroup in get_naba_ML', inode
          if(nn.lt.prim%groups_on_node) then
             amat(nn+1)%offset=amat(nn)%offset+amat(nn)%part_nabs
             amat(nn+1)%nd_offset=amat(nn)%nd_offset+amat(nn)%part_nd_nabs
          endif
       else
          amat(nn)%n_atoms = 0
          if(nn.lt.prim%groups_on_node) then
             amat(nn+1)%offset=amat(nn)%offset+amat(nn)%part_nabs
             amat(nn+1)%nd_offset=amat(nn)%nd_offset+amat(nn)%part_nd_nabs
          endif
          amat(nn)%i_acc(1) = 1
          amat(nn)%i_nd_acc(1) = 1
       endif ! End if(prim%nm_nodgroup>0)
    enddo ! End part_on_node

    ! Now create length of the matrix
    amat%length=0
    do nn=1,prim%groups_on_node
       amat(1)%length=amat(1)%length+amat(nn)%part_nd_nabs
    enddo
    amat(2:prim%groups_on_node)%length = amat(1)%length
    return
  end subroutine get_naba_ML
!!***


!!****f* mlff_module/get_triplet_ML *
!! modified from matrix_elements_module/get_naba
!!
!!  NAME
!!   get_triplet_ML
!!  USAGE
!!
!!  PURPOSE
!!   Creates triplet lists for a matrix range rcut
!!  INPUTS
!!    integer(integ) :: part_on_node
!!    real(double) :: rcut
!!    type(matrix_ML), dimension (:) :: amat
!!  USES
!!
!!  AUTHOR
!!   JianBo.Lin
!!  CREATION DATE
!!   2022/07/20
!!  MODIFICATION HISTORY
!!   2022/07/27 J.Lin
!!    Modified to get information of triplets for machine learing
!!  SOURCE
!!
  subroutine get_triplet_ML(amat,part_on_node,rcut)

    ! Module usage
    use datatypes

    use basic_types
    use units
    use numbers,   ONLY: pi,very_small
    use GenComms, ONLY: cq_abort, inode, ionode, my_barrier
    use mlff_type, ONLY: matrix_ML

    implicit none

    ! Passed variables
    type(matrix_ML), dimension (:), intent(inout) :: amat
    integer(integ), intent(in) :: part_on_node
    real(double), intent(in) :: rcut

    ! Local variables
    integer(integ) :: inp,nn
    integer(integ) :: i, j, k, ist_j, ist_k, ist_ijk
    real(double) :: xij, yij, zij, xik, yik, zik, rjk, rjk_2,rcut_2

    ! After we have neighbor information and descriptor information
    ! Get information of triplets
    ! In this subrputine, Units: Length=Bohr
    rcut_2=rcut * rcut
    call my_barrier()

    ! loop over all atom pairs (atoms in primary set, max. cover set) -
    inp=1  ! Indexes primary atoms
    do nn=1,part_on_node ! Partitions in primary set
      if (amat(nn)%n_atoms .gt. 0) then ! Are there atoms ?
      ! collect accumulation of three body terms, triplets
      amat(nn)%triplet_acc(1)=1
      do i=1,amat(nn)%n_atoms  ! Loop over atoms in partition

        do j=1,  amat(nn)%n_nab(i)
            ist_j = amat(nn)%i_acc(i)+j-1
            xij=amat(nn)%dx(ist_j)
            yij=amat(nn)%dy(ist_j)
            zij=amat(nn)%dz(ist_j)

          amat(nn)%n_triplet(ist_j)=0
          do k=j+1,  amat(nn)%n_nab(i)
            ist_k = amat(nn)%i_acc(i)+k-1
            xik=amat(nn)%dx(ist_k)
            yik=amat(nn)%dy(ist_k)
            zik=amat(nn)%dz(ist_k)
            rjk_2 = (xij-xik)**2+(yij-yik)**2+(zij-zik)**2

            if (rjk_2 >very_small .and. rjk_2 < rcut_2) then
              amat(nn)%n_triplet(ist_j)=amat(nn)%n_triplet(ist_j)+1
              ist_ijk=amat(nn)%triplet_acc(ist_j)+amat(nn)%n_triplet(ist_j)-1
              amat(nn)%ist_triplet(ist_ijk)=ist_k
            end if ! rjk < rcut
          end do ! k, three-body
          if(ist_j.lt.amat(nn)%part_nabs) then
             amat(nn)%triplet_acc(ist_j+1)=amat(nn)%triplet_acc(ist_j)+amat(nn)%n_triplet(ist_j)
          endif
          amat(nn)%part_triplets = amat(nn)%part_triplets+amat(nn)%n_triplet(ist_j)
          end do ! j, two-body
        inp=inp+1  ! Indexes primary-set atoms
        enddo ! End i in prim%nm_nodgroup

      else
        ! if no atoms, do nothing, and turn to next partition
        write(*, *) 'Warning: no atom in partition,', nn, amat(nn)%n_atoms
      end if
    enddo ! End part_on_node
  end subroutine get_triplet_ML
!!***


!!****f* mlff_type_module/get_feature_acsf2b *
!!
!!  NAME
!!   get_feature_acsf2b
!!  USAGE
!!
!!  PURPOSE
!!   calculate ml force and energy for each atom
!!  INPUTS
!!   type(primary_set) :: prim
!!   type(cover_set) :: gcs
!!   type(matrix), dimension(:) :: amat
!!   real(double) :: rcut
!!   type(acsf_param), dimension(:) :: descriptor_params
!!  USES
!!
!!  AUTHOR
!!   JianBo.Lin
!!  CREATION DATE
!!   2022/07/20
!!  MODIFICATION HISTORY
!!   2022/07/27 JL
!!    Modified to get pair information for machine learing
!!  SOURCE
!!
  subroutine calculate_EandF_acsf2b(prim,amat,amat_features_ML,descriptor_params)

    ! Module usage
    use datatypes
    use basic_types
    use matrix_module
    use GenComms, ONLY: cq_abort, inode, ionode, my_barrier
    use global_module, ONLY: id_glob, species_glob, sf, nlpf, paof, napf, numprocs,&
            flag_stress, flag_full_stress
    use group_module, ONLY: parts
    use species_module, ONLY: nsf_species, nlpf_species, npao_species, napf_species

    use energy, ONLY: ml_energy
    use mlff_type, ONLY:matrix_ML,features_ML,acsf_param

    implicit none

    ! Passed variables
    type(primary_set), intent(in) :: prim
    type(matrix_ML), dimension (:), intent(in) :: amat
    type(features_ML), dimension(:), intent(inout) :: amat_features_ML
    type(acsf_param), dimension(:), intent(in) :: descriptor_params

    ! Local variables
    integer :: inp, ip_process
    integer :: nn,i,np, i_species,ia_glob
    real(double) :: rx,ry,rz

    ! loop over all atom pairs (atoms in primary set, max. cover set) -
    inp=1  ! Indexes primary atoms
    do nn=1,prim%groups_on_node ! Partitions in primary set
       !pair check!write(*,*) 'nn loop, inode=',inode
       if(prim%nm_nodgroup(nn).gt.0) then  ! Are there atoms ?
          do i=1,prim%nm_nodgroup(nn)  ! Loop over atoms in partition
             ia_glob=prim%ig_prim(prim%nm_nodbeg(nn)+i-1)
             i_species = amat(nn)%i_species(i)
             ! TODO: fpx(i, :) -> fpx(:, i)
             ml_force(1, ia_glob) = dot_product(amat_features_ML(nn)%fpx(:, i), descriptor_params(i_species)%coef)
             ml_force(2, ia_glob) = dot_product(amat_features_ML(nn)%fpy(:, i), descriptor_params(i_species)%coef)
             ml_force(3, ia_glob) = dot_product(amat_features_ML(nn)%fpz(:, i), descriptor_params(i_species)%coef)

             ! Get ml_stress for each NB atom
             if(flag_stress) then
                 rx=prim%xprim(inp)
                 ry=prim%yprim(inp)
                 rz=prim%zprim(inp)
                 if(flag_full_stress) then
                    ! vec_Fi*vec_ri
                     ml_stress(1,1) = ml_stress(1,1) + ml_force(1, ia_glob) * rx
                     ml_stress(1,2) = ml_stress(1,2) + ml_force(1, ia_glob) * ry
                     ml_stress(1,3) = ml_stress(1,3) + ml_force(1, ia_glob) * rz

                     ml_stress(2,1) = ml_stress(2,1) + ml_force(2, ia_glob) * rx
                     ml_stress(2,2) = ml_stress(2,2) + ml_force(2, ia_glob) * ry
                     ml_stress(2,3) = ml_stress(2,3) + ml_force(2, ia_glob) * rz

                     ml_stress(3,1) = ml_stress(3,1) + ml_force(3, ia_glob) * rx
                     ml_stress(3,2) = ml_stress(3,2) + ml_force(3, ia_glob) * ry
                     ml_stress(3,3) = ml_stress(3,3) + ml_force(3, ia_glob) * rz
                 else
                     ml_stress(1,1) = ml_stress(1,1) + ml_force(1, ia_glob) * rx
                     ml_stress(2,2) = ml_stress(2,2) + ml_force(2, ia_glob) * ry
                     ml_stress(3,3) = ml_stress(3,3) + ml_force(3, ia_glob) * rz
                 end if
             end if
             inp=inp+1  ! Indexes primary-set atoms
          enddo ! End prim%nm_nodgroup
          !pair check!write(*,*) 'after prim%nm_nodgroup', inode
       else
            write(*, *) 'Warning: No atoms in this partition', inode, nn
       endif ! End if(prim%nm_nodgroup>0)
    enddo ! End part_on_node
  end subroutine calculate_EandF_acsf2b
!!***


!!****f* mlff_type_module/calculate_EandF_split2b3b *
!!
!!  NAME
!!   calculate_EandF_split2b3b
!!  USAGE
!!
!!  PURPOSE
!!   calculate ml force and energy for each atom
!!  INPUTS
!!   type(primary_set) :: prim
!!   type(cover_set) :: gcs
!!   type(matrix), dimension(:) :: amat
!!   real(double) :: rcut
!!   type(acsf_param), dimension(:) :: descriptor_params
!!  USES
!!
!!  AUTHOR
!!   JianBo.Lin
!!  CREATION DATE
!!   2022/07/20
!!  MODIFICATION HISTORY
!!   2022/07/27 JL
!!    Modified to get pair information for machine learing
!!  SOURCE
!!
  subroutine calculate_EandF_split2b3b(prim,amat,amat_features_ML,descriptor_params)

    ! Module usage
    use datatypes
    use basic_types
    use matrix_module
    use GenComms, ONLY: cq_abort, inode, ionode, my_barrier
    use global_module, ONLY: id_glob, species_glob, sf, nlpf, paof, napf, numprocs,&
            flag_stress, flag_full_stress
    use group_module, ONLY: parts
    use species_module, ONLY: nsf_species, nlpf_species, npao_species, napf_species

    use energy, ONLY: ml_energy
    use mlff_type, ONLY:matrix_ML,features_ML,split_param

    implicit none

    ! Passed variables
    type(primary_set), intent(in) :: prim
    type(matrix_ML), dimension (:), intent(in) :: amat
    type(features_ML), dimension(:), intent(inout) :: amat_features_ML
    type(split_param), dimension(:), intent(in) :: descriptor_params

    ! Local variables
    integer :: inp, ip_process
    integer :: nn,i,np, i_species,ia_glob
    real(double) :: rx,ry,rz

    ! loop over all atom pairs (atoms in primary set, max. cover set) -
    inp=1  ! Indexes primary atoms
    do nn=1,prim%groups_on_node ! Partitions in primary set
       !pair check!write(*,*) 'nn loop, inode=',inode
       if(prim%nm_nodgroup(nn).gt.0) then  ! Are there atoms ?
          do i=1,prim%nm_nodgroup(nn)  ! Loop over atoms in partition
             ia_glob=prim%ig_prim(prim%nm_nodbeg(nn)+i-1)
             i_species = amat(nn)%i_species(i)

             ml_force(1, ia_glob) = dot_product(amat_features_ML(nn)%fpx(:, i), descriptor_params(i_species)%coef)
             ml_force(2, ia_glob) = dot_product(amat_features_ML(nn)%fpy(:, i), descriptor_params(i_species)%coef)
             ml_force(3, ia_glob) = dot_product(amat_features_ML(nn)%fpz(:, i), descriptor_params(i_species)%coef)

             ! Get ml_stress for each NB atom
             if(flag_stress) then
                 rx=prim%xprim(inp)
                 ry=prim%yprim(inp)
                 rz=prim%zprim(inp)
                 if(flag_full_stress) then
                    ! vec_Fi*vec_ri
                     ml_stress(1,1) = ml_stress(1,1) + ml_force(1, ia_glob) * rx
                     ml_stress(1,2) = ml_stress(1,2) + ml_force(1, ia_glob) * ry
                     ml_stress(1,3) = ml_stress(1,3) + ml_force(1, ia_glob) * rz

                     ml_stress(2,1) = ml_stress(2,1) + ml_force(2, ia_glob) * rx
                     ml_stress(2,2) = ml_stress(2,2) + ml_force(2, ia_glob) * ry
                     ml_stress(2,3) = ml_stress(2,3) + ml_force(2, ia_glob) * rz

                     ml_stress(3,1) = ml_stress(3,1) + ml_force(3, ia_glob) * rx
                     ml_stress(3,2) = ml_stress(3,2) + ml_force(3, ia_glob) * ry
                     ml_stress(3,3) = ml_stress(3,3) + ml_force(3, ia_glob) * rz
                 else
                     ml_stress(1,1) = ml_stress(1,1) + ml_force(1, ia_glob) * rx
                     ml_stress(2,2) = ml_stress(2,2) + ml_force(2, ia_glob) * ry
                     ml_stress(3,3) = ml_stress(3,3) + ml_force(3, ia_glob) * rz
                 end if
             end if
             inp=inp+1  ! Indexes primary-set atoms
          enddo ! End prim%nm_nodgroup
          !pair check!write(*,*) 'after prim%nm_nodgroup', inode
       else
            write(*, *) 'Warning: No atoms in this partition', inode, nn
       endif ! End if(prim%nm_nodgroup>0)
    enddo ! End part_on_node
  end subroutine calculate_EandF_split2b3b
!!***

!!****f* mlff_type_module/get_force_ML *
!! modified from matrix_elements_module/get_naba
!!
!!  NAME
!!   get_force_ML
!!  USAGE
!!
!!  PURPOSE
!!   Creates neighbour lists for a matrix range rcut
!!  INPUTS
!!   type(primary_set) :: prim
!!   type(cover_set) :: gcs
!!   type(matrix), dimension(:) :: amat
!!   real(double) :: rcut
!!  USES
!!
!!  AUTHOR
!!   JianBo.Lin
!!  CREATION DATE
!!   2022/07/20
!!  MODIFICATION HISTORY
!!   2022/07/27 JL
!!    Modified to get pair information for machine learing
!!  SOURCE
!!
  subroutine get_EandF_ML(prim,gcs,amat_ML,amat_features_ML,rcut)

    ! Module usage
    use datatypes
    use basic_types
    use matrix_module
    use GenComms, ONLY: cq_abort, inode, ionode, my_barrier
    use global_module, ONLY: id_glob, species_glob, sf, nlpf, paof, napf, numprocs
    use group_module, ONLY: parts
    use species_module, ONLY: nsf_species, nlpf_species, npao_species, napf_species

    use mlff_type

    implicit none

    ! Passed variables
    type(primary_set), intent(in) :: prim
    type(cover_set), intent(in) :: gcs
    real(double), intent(in) :: rcut
    type(matrix_ML), dimension (:), allocatable, intent(inout) :: amat_ML
    type(features_ML), dimension (:), allocatable, intent(inout) :: amat_features_ML

    ! Local variables
    integer                    :: file_id=2023, stat
    character(len=80)          :: filename
    character(len=80)          :: flag_descriptortype, comment

    integer :: inp, cumu_ndims, neigh_spec, ip_process,feature_dim
    integer :: nn,j,np,ni,ist
    real(double) :: rcutsq,dx,dy,dz,rij2
    real(double), parameter :: tol=1.0e-8_double

    ! Todo: filename
    filename='ml_pot.txt'
    ! After we have neighbor information and descriptor information
    ! read and set
    open(unit=file_id,file=filename,iostat=stat,status='old')
    read(file_id,*) comment
    read(file_id,*) flag_descriptortype
    close(file_id)

    ! Todo: check flag_descriptortype
    if (inode == ionode) then
        write(*,*) 'flag_descriptortype:', flag_descriptortype
    end if


    if (params_ML%descriptor_type == 'split2b3b' .or. params_ML%descriptor_type == 'Split2b3b') then
        !call read_split(filename, descriptor_params_split)
        feature_dim = params_ML%split(1)%dim_coef
        call allocate_features_ML(amat_ML, amat_features_ML, prim%groups_on_node, feature_dim)
        call get_feature_split(prim,gcs,amat_ML,amat_features_ML,rcut,params_ML%split)
        call calculate_EandF_split2b3b(prim,amat_ML,amat_features_ML,params_ML%split)

        !call deallocate_split_param(descriptor_params_split)
        !call clean_up_ml_split(amat_ML,amat_features_ML,descriptor_params_split)
    elseif (params_ML%descriptor_type == 'acsf2b' .or. params_ML%descriptor_type == 'ACSF2b') then
        !call read_acsf2b(filename, descriptor_params_acsf2b)
        feature_dim = params_ML%acsf(1)%dim_coef
        call allocate_features_ML(amat_ML, amat_features_ML, prim%groups_on_node, feature_dim)
        call get_feature_acsf2b(prim,gcs,amat_ML,amat_features_ML,rcut,params_ML%acsf)
        call calculate_EandF_acsf2b(prim,amat_ML,amat_features_ML,params_ML%acsf)

        !call deallocate_acsf2b_param(descriptor_params_acsf2b)
        !call clean_up_ml_acsf(amat_ML,amat_features_ML,descriptor_params_acsf2b)
    else
        call cq_abort('flag_descriptortype: (Error) no such descriptor type at present')
    end if

    !call deallocate_features_ML(amat_ML,amat_features_ML,prim%groups_on_node)
    !call deallocate_matrix_ML(amat_ML,prim%groups_on_node,parts%mx_mem_grp)
  end subroutine get_EandF_ML
!!***


!!****f* mlff_module/get_MLFF *
!!
!!  NAME
!!   get_MLFF
!!  PURPOSE
!!   Calculate machine learning forces and energies
!!  AUTHOR
!!   JianBo.Lin
!!  CREATION DATE
!!   2022/07/27
!!  MODIFICATION HISTORY
!!
!!  SOURCE
!!
  subroutine get_MLFF

    use datatypes
    use units
    use numbers
    use dimens, ONLY: atomicnum, r_ML_des, r_super_x, r_super_y, r_super_z
    use cover_module, ONLY: make_cs, ML_CS, make_iprim, send_ncover
    use GenComms, ONLY: cq_abort, inode, ionode, gsum, my_barrier
    use group_module, ONLY: parts
    use global_module, ONLY: io_lun, numprocs, ni_in_cell, &
            species_glob, id_glob, flag_stress, flag_full_stress,&
            x_atom_cell, y_atom_cell, z_atom_cell
    use primary_module, ONLY: bundle
    use io_module, ONLY: get_file_name
    use species_module, ONLY: n_species, species_label
    use ion_electrostatic, ONLY: partition_distance
    use mlff_type
    use matrix_module, ONLY: matrix, deallocate_matrix
    use force_module, ONLY: tot_force,stress
    use mpi

    implicit none

    ! local variables
    integer :: i, j, ip, ia, np, nj, ia_glob, j_glob,nn, ist
    integer :: stat !, Z_i, Z_j
    !integer :: nppx, nppy, nppz, ind_part, nccx, nccy, nccz, m1, m2, &
    !     m3, m, ig_atom_beg
    integer, allocatable :: neighbour_part(:)
    real(double) :: real_cell_vec(3, 3), part_cell_vec(3, 3), &
            part_cell_dual(3, 3)
    real(double) :: t1,t2
    type(matrix), dimension (:), allocatable :: amat
    type(matrix_ML), dimension (:), allocatable :: amat_ML
    type(features_ML), dimension (:), allocatable :: amat_features_ML

    if (inode== ionode) then
        write(*,*) 'We are in get_MLFF'
    end if

    t1=MPI_wtime()
    !! Make covering set
    call make_cs(inode-1, r_ML_des, ML_CS, parts, bundle, ni_in_cell, &
                   x_atom_cell, y_atom_cell, z_atom_cell)
    !call set_para_ML

    ! Initialisation of the machine learning energy and the atomic forces
    allocate ( ml_force(3, ni_in_cell), STAT = stat )
    if (stat .NE. 0) call cq_abort("Error allocating ml_force: ", ni_in_cell)

    ml_energy = 0.0_double
    ml_force  = 0.0_double
    if(flag_stress) ml_stress = zero

    ! Make ML cover set and send
    call make_iprim(ML_CS, bundle)
    call send_ncover(ML_CS, inode) ! comm
    call my_barrier

    allocate ( neighbour_part(ML_CS%ng_cover), STAT=stat )
    if (stat .NE. 0) call cq_abort("Eror allocating neighbour_list in get_MLFF: ", ML_CS%ng_cover)
    allocate ( amat(bundle%groups_on_node), STAT=stat)
    if (stat .NE. 0) call cq_abort("Eror allocating amat in get_MLFF: ", bundle%groups_on_node)
    allocate ( amat_ML(bundle%groups_on_node), STAT=stat)
    if (stat .NE. 0) call cq_abort("Eror allocating amat_ML in get_MLFF: ", bundle%groups_on_node)
    allocate ( amat_features_ML(bundle%groups_on_node), STAT=stat)
    if (stat .NE. 0) call cq_abort("Eror allocating amat_features_ML in get_MLFF: ", bundle%groups_on_node)

    real_cell_vec(1, 1) = r_super_x
    real_cell_vec(1, 2) = zero
    real_cell_vec(1, 3) = zero
    real_cell_vec(2, 1) = zero
    real_cell_vec(2, 2) = r_super_y
    real_cell_vec(2, 3) = zero
    real_cell_vec(3, 1) = zero
    real_cell_vec(3, 2) = zero
    real_cell_vec(3, 3) = r_super_z

    ! primitive translation vectors for partitions
    do i = 1, 3
      part_cell_vec(1, i) = real_cell_vec(1, i) / real(parts%ngcellx, double)
      part_cell_vec(2, i) = real_cell_vec(2, i) / real(parts%ngcelly, double)
      part_cell_vec(3, i) = real_cell_vec(3, i) / real(parts%ngcellz, double)
    enddo

    ! Initialise the partition distance routine
    call partition_distance(.true., part_cell_vec, part_cell_dual)
    !! TODO:If necessary my_barrier
    !call my_barrier()

    ! Kernel of the ML scheme
    if (inode== ionode) then
        write(*,*) 'before get_naba_ML:: bundle%n_prim=', bundle%n_prim, bundle%groups_on_node
    end if

    ! Initialise the ML matrix and Covering set
    if(bundle%n_prim.gt.0) then
        call matrix_ini_ML(parts,bundle,ML_CS,amat,amat_ML, r_ML_des,inode)
    end if
    !! TODO:If necessary my_barrier
    !call my_barrier()

    !! Check after get naba ML
    if (inode == ionode) then
        write(*,*) '################# After get_naba_ML #################'
        write(*,*) 'sign,  inode,   id_partition, iatom_species'
    end if
    !! Check ia glob id
    do nn=1, bundle%groups_on_node
        if(bundle%nm_nodgroup(nn).gt.0) then
            do i=1, amat_ML(nn)%n_atoms
                ia_glob=bundle%ig_prim(bundle%nm_nodbeg(nn)+i-1)
                !check!write(*, *) 'i_species',inode, nn, 'ia_glob1', ia_glob,amat_ML(nn)%i_species(i)
            end do !part_nabs
        end if
    end do

    !check get_naba_ML!if (inode == ionode) then
        !check get_naba_ML!write(*,101) 'inode:','i_par:','i_id:','i_glob2:', 'j_par:','j_id:','j_glob2:', &
        !check get_naba_ML!    'j_species', 'r(Angstrom)'
    !check get_naba_ML!end if

    !! Check nb list
    if (inode == ionode) then
        do nn=1, bundle%groups_on_node
            if(bundle%nm_nodgroup(nn).gt.0) then
                do i=1, amat_ML(nn)%n_atoms
                    ia_glob=bundle%ig_prim(bundle%nm_nodbeg(nn)+i-1)

                    !check!if (ia_glob==1) then
                        !check!write(1984,*) 'this is checking first atom in x direction: before'
                        !check!write(1984,*) amat_features_ML(nn)%fpx(i, :)
                        !check!write(1984,*) 'fx= ',ml_force(1, ia_glob)
                    !check!end if

                    do j=1, amat_ML(nn)%n_nab(i)
                        ist=amat_ML(nn)%i_acc(i)+j-1

                        np=amat_ML(nn)%i_part(ist)
                        nj=amat_ML(nn)%i_seq(ist)
                        j_glob =  id_glob( parts%icell_beg(ML_CS%lab_cell(np)) +nj-1 )

                        !write(*, 100) inode, nn,i, ia_glob,np,nj, j_glob,&
                        !        amat_ML(nn)%j_species(ist),amat_ML(nn)%radius(ist) * BohrToAng
                        !check!write(*, 102) inode, nn,i, ia_glob,np,nj, j_glob,&
                                !check!amat_ML(nn)%j_species(ist),amat_ML(nn)%radius(ist) * BohrToAng
                    end do
                end do !part_nabs
            else
                write(*, *) 'Warning: No atoms in this partition', inode, nn
            end if

        end do
        write(*,*) '################# END Check get_naba_ML #################'
    end if
    !! End of Check after get naba ML
    t2=MPI_wtime()
    if (inode==ionode) &
      write(io_lun,*) 'Time before get_E_and_F_ML in get_MLFF:', t2-t1

    !! mpi_time()
    !! get_EandF_ML: read infomation of ML,
    t1=MPI_wtime()
    call get_EandF_ML(bundle,ML_CS,amat_ML,amat_features_ML,r_ML_des)
    ! remove in future
    call my_barrier()
    t2=MPI_wtime()
    if (inode==ionode) &
      write(io_lun,*) 'Time get_E_and_F_ML in get_MLFF:', t2-t1


    !! mpi_time()
    t1=MPI_wtime()
    call gsum(ml_energy)
    call gsum(ml_force, 3, ni_in_cell)
    tot_force = ml_force * BohrToAng / HaToeV
    if(flag_stress) then
       call gsum(ml_stress,3,3)
       stress = -ml_stress / HaToeV
    end if
    t2=MPI_wtime()
    if (inode==ionode) &
      write(io_lun,*) 'Time after get_E_and_F_ML in get_MLFF:', t2-t1

    !! check force and feature info after operate mlff
    if (inode == ionode) then
        do nn=1, bundle%groups_on_node
            if(bundle%nm_nodgroup(nn).gt.0) then
                do i=1, amat_ML(nn)%n_atoms
                    ia_glob=bundle%ig_prim(bundle%nm_nodbeg(nn)+i-1)

                    !!if (ia_glob==1 .or. ia_glob==71) then
                    if (ia_glob==1) then
                        write(1984,*) 'this is checking first atom in x direction: after, ia_glob:',ia_glob
                        do j=1, 40
                            write(1984,*) 'No' !,amat_features_ML(nn)%fpx(i, j)
                        end do
                        write(1984,*) 'fx= ',ml_force(1, ia_glob)
                    end if
                end do !i
            else
                write(*, *) 'Warning: No atoms in this partition', inode, nn
            end if ! no atoms in partition?
        end do !nn

        write(1988,*) ml_energy
        do i=1, ni_in_cell
            write(1987,103) i, ml_force(:, i)
        end do
    end if

    ! Deallocate feature and matrix
    call deallocate_features_ML(amat_ML,amat_features_ML,bundle%groups_on_node)
    call deallocate_matrix_ML(amat_ML,bundle%groups_on_node,parts%mx_mem_grp)
    !deallocate(amat)
    !call deallocate_matrix(amat,bundle%groups_on_node)

    deallocate (neighbour_part, STAT=stat )
    if (stat .NE. 0) call cq_abort("Eror deallocating neighbour_part: ", stat)
    deallocate(ml_force)

    !!call end_comms()
    if (inode == ionode) then
        write(*,*) '################# END Check get_MLFF #################'
    end if


    100 format('inode: ',i5,' i_par: ',i5,' i_id: ',i5,' i_glob2: ', i5,' j_par: ',i5,' j_id: ',i5,' j_glob2: ', i5, &
            ' j_species: ', i5,' radius(Angstrom): ',f12.8)
    101 format(a8,a8,a8,a8, a8,a8,a8, a10,a12)
    102 format(i8,i8,i8,i8, i8,i8,i8, i10,f12.8)
    103 format(i8,3f16.8)

    !stop
    !!call gsum(disp_energy)
    !!call gsum(disp_force, 3, ni_in_cell)
    !!disp_energy = - disp_energy * s_6 * half ! CANCEL THE DOUBLE-COUNTING
    !!disp_force = s_6 * disp_force
    !!if(flag_stress) then
    !!   call gsum(ml_stress,3,3)
    !!   ml_stress = -half * s_6 * ml_stress  ! Double counting again
    !!end if

    ! Works only when you want to get the dispersion
    !!if (flag_only_dispersion) then
    !!  if (inode == ionode) then
    !!    write (io_lun, '(/a)') "***** Dispersion Information: DFT-D2 *****"
    !!    write (io_lun,16) en_conv*disp_energy, en_units(energy_units)
    !!    write (io_lun,fmt='(/,20x,"Disp. Force on atoms (",a2,"/",a2,")"/)') en_units(energy_units), d_units(dist_units)
    !!    write (io_lun,fmt='("Atom        X               Y               Z")')
    !!    do i = 1, ni_in_cell
    !!        write (io_lun, 109) i, (for_conv*disp_force(j,i), j=1,3)
    !!     enddo
    !!     write(io_lun,'(/8x,"Dispersion stress: ",3f15.8)') ml_stress
    !!  endif
    !!endif

    110 format('E_ML:',f25.15,' ',a2)
    111 format('Force on atom ', i9)
    112 format(i10, 3f15.10)
    return
  end subroutine get_MLFF
!!***

end module mlff
