!!****h* Conquest/store_matrix_module
!!  NAME
!!   store_matrix_module
!!  PURPOSE
!!    Stores the data and information of a matrix for a set of atomic positions.
!!   to use them for a different set of atomic positions. (ex. during the MD)
!!  USES
!!
!!  AUTHOR
!!   Tsuyoshi Miyazaki
!!  CREATION DATE
!!   2016/09/29
!!  MODIFCATION
!!   - Added derived data types for XL-BOMD
!!
module store_matrix

  use datatypes
  use global_module, ONLY: numprocs
  use input_module, ONLY: io_assign, io_close
  use GenComms, ONLY: inode, ionode, cq_abort

  implicit none

  type matrix_store
    character(len=80) :: name
    !integer :: inode  ! not necessary 
    integer :: n_prim
    integer :: matrix_size
   ! Size : n_prim
    integer, pointer :: nsf_spec_i(:)
    integer, pointer :: idglob_i(:)
    integer, pointer :: jmax_i(:)
    integer, pointer :: jbeta_max_i(:)
    integer, pointer :: ibeg_data_matrix(:)  ! initial address of data_matrix for (iprim)th atom
   ! Size : (n_prim, jmax_i_max)
    integer, pointer :: idglob_j(:,:)
    integer, pointer :: beta_j(:,:)
   ! Size : (3, jmax_i_max, iprim)
    real(double), pointer :: vec_Rij(:,:,:)
   ! Size : matrix_size
    real(double), pointer :: data_matrix(:)
  end type matrix_store

  character(80),private :: RCSid = "$Id$"

contains
  ! -----------------------------------------------------------------------
  ! Subroutine dump_matrix2
  ! -----------------------------------------------------------------------

  !!****f* io_module2/dump_matrix2 *
  !!
  !!  NAME
  !!    dump_matrix2
  !!  USAGE
  !!
  !!  PURPOSE
  !!    Stores the information as for the primary atoms and its
  !!    neighbours first, then the matrix elements are followed.
  !!    Note that matrix data are the ones obtained at the previous (old)
  !!    step - e.g. L_{i.OLD,alpha | j.OLD,beta}.
  !!
  !!    imax:           the number of primary atoms "i" in the current node "inode".
  !!    alpha_i():      the number of orbitals for each primary set of atom "i"
  !!                    in the current node "inode".
  !!    jmax():         the number of neighbour-atoms "j" for each atom "i" above.
  !!    jbeta_max():    the number of (jmax x beta) for each atom "i".
  !!    beta_j():       the number of orbitals for each neighbour-atom "j".
  !!    idglob_j():     the global labelling of neighbour-atoms "j".
  !!
  !!  INPUTS
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   Michiaki Arita
  !!  CREATION DATE
  !!   2013/08/21
  !!  MODIFICATION
  !!   2016/10/03 Tsuyoshi Miyazaki
  !!      1. Changed to use the derivative type (matrix_store) for tidying up MD continuation
  !!      2. Moved from io_module2
  !!  SOURCE
  !!

   subroutine dump_matrix2(stub,matA,range)
    use GenComms, ONLY: inode, ionode, cq_abort
    use global_module, ONLY: numprocs, id_glob
    use io_module, ONLY: get_file_name
    implicit none
    character(len=*),intent(in) :: stub
    integer,intent(in) :: matA
    integer,intent(in) :: range
    type(matrix_store):: tmp_matrix_store

    integer :: lun, iprim, nprim, jmax, jj, ibeg, jbeta_alpha, len
    character(20) :: file_name

    ! set_matrix_store : build tmp_matrix_store
    call set_matrix_store(stub,matA,range,tmp_matrix_store)

    ! Actual Dump (from dump_matrix2)
     ! First, make a file based upon the node ID.
     call get_file_name(stub//'matrix2',numprocs,inode,file_name)
     call io_assign(lun)
    open (lun,file=file_name)

      ! 1. node ID, no. of PS of atoms "i".
        nprim=tmp_matrix_store%n_prim
        write (lun,*) inode, nprim
      ! 2. no. of alpha for each "i".
        write (lun,*) tmp_matrix_store%nsf_spec_i(1:nprim)
      ! NOTE: The followings are written out with neighbour-labelling
      ! 3. no. of the neighbours "j" for each "i".
      ! 4. no. of "neighbour-j x beta" for each "i".
        write (lun,*) tmp_matrix_store%jmax_i(1:nprim)
        write (lun,*) tmp_matrix_store%jbeta_max_i(1:nprim)

     !I will change the order of dumping in the following, later.   2016/09/30: TM@UCL
      if(nprim .GT. 1) then
       do iprim=1,nprim
          jmax = tmp_matrix_store%jmax_i(iprim)
          write (lun,*) tmp_matrix_store%idglob_i(iprim)
          write (lun,*) tmp_matrix_store%beta_j(1:jmax,iprim)
          write (lun,*) tmp_matrix_store%idglob_j(1:jmax,iprim)
         do jj=1,jmax
          write (lun,*) tmp_matrix_store%vec_Rij(1:3,jj,iprim)
         enddo !jj=1,jmax
         if(iprim < nprim) then
           len = tmp_matrix_store%ibeg_data_matrix(iprim+1)-tmp_matrix_store%ibeg_data_matrix(iprim)
         else
           len = tmp_matrix_store%matrix_size-tmp_matrix_store%ibeg_data_matrix(iprim)+1
         endif
           ibeg = tmp_matrix_store%ibeg_data_matrix(iprim)
         do jbeta_alpha = 1, len
          write (lun,fmt='(f25.18)') tmp_matrix_store%data_matrix(ibeg+jbeta_alpha-1)
         enddo !jbeta_alpha = 1, len
       enddo !iprim=1,nprim
      endif  ! (nprim .GT. 1) 

    ! Close the file in the end.
    call io_close(lun)
    
    ! free_matrix_store : free tmp_matrix_store
    call free_matrix_store(tmp_matrix_store)

    return
   end subroutine dump_matrix2

  ! -----------------------------------------------------------------------
  ! Subroutine set_matrix_store  & free_matrix_store
  ! -----------------------------------------------------------------------
  !!****f* store_matrix_module/set_matrix_store *
  !!
  !!  NAME
  !!    set_matrix_store
  !!  USAGE
  !!
  !!  PURPOSE
  !!    Stores the information as for the primary atoms and its
  !!    neighbours first, then the matrix elements are followed.
  !!    Made from dump_matrix2 written by M. Arita
  !!
  !!  INPUTS
  !!
  !!  USES
  !!
  !!  AUTHOR
  !!   Tsuyoshi Miyazaki/Michiaki Arita
  !!  CREATION DATE
  !!   2016/10/03
  !!  MODIFICATION
  !!
  !!  SOURCE
  !!

   subroutine set_matrix_store(stub,matA,range,matinfo)
    use numbers
    use global_module, ONLY: id_glob, species_glob
    use species_module, ONLY: nsf_species
    use group_module, ONLY: parts
    use primary_module, ONLY: bundle
    use cover_module, ONLY: BCS_parts
    use matrix_data, ONLY: mat
    use mult_module, ONLY: mat_p

    implicit none
    character(len=*),intent(in) :: stub
    integer,intent(in) :: matA
    integer,intent(in) :: range
    type(matrix_store), intent(out) :: matinfo
    !Local
    integer :: nprim, istat
    integer :: np, ni, iprim, atom_num, jmax_i_max
    integer :: ist, n_beta, neigh, gcspart, ibeg, j_global_part, len

    !name & n_prim
      matinfo%name=stub
      matinfo%n_prim=bundle%n_prim
    !allocation 1 
      nprim=matinfo%n_prim   ! local variable
      allocate(matinfo%nsf_spec_i(nprim), matinfo%idglob_i(nprim), matinfo%jmax_i(nprim), &
               matinfo%jbeta_max_i(nprim), matinfo%ibeg_data_matrix(nprim), STAT=istat)
      if(istat .NE. 0) call cq_abort('Allocation 1 in set_matrix_store',istat,nprim)
    
      matinfo%nsf_spec_i(1:nprim) =nsf_species(bundle%species(1:nprim))

   !!set up : jmax_i, jbeta_max_i

    iprim = 0
    jmax_i_max = -1
    matinfo%jbeta_max_i(:) = 0

    do np = 1, bundle%groups_on_node
      if (bundle%nm_nodgroup(np).GT.0) then
        do ni = 1, bundle%nm_nodgroup(np)
           iprim = iprim + 1
           atom_num = bundle%nm_nodbeg(np)+ni-1
           if(iprim .NE. atom_num) call cq_abort('iprim .ne. atom_num ?',iprim,atom_num)
           matinfo%idglob_i(iprim) = bundle%ig_prim(iprim)
           matinfo%jmax_i(iprim) =  mat(np,range)%n_nab(ni)
           if (matinfo%jmax_i(iprim).GT.jmax_i_max) jmax_i_max = matinfo%jmax_i(iprim)
          do neigh = 1, matinfo%jmax_i(iprim)
            ist = mat(np,range)%i_acc(ni) + neigh - 1 ! Neighbour-labeling
            n_beta = mat(np,range)%ndimj(ist)
            matinfo%jbeta_max_i(iprim) = matinfo%jbeta_max_i(iprim) + n_beta
          enddo !(neigh)
        enddo !(ni)
      endif
    enddo !(np)

    !allocation 2 
      allocate(matinfo%idglob_j(jmax_i_max,nprim), matinfo%beta_j(jmax_i_max,nprim), STAT=istat)
      if(istat .NE. 0) call cq_abort('Allocation 2 in set_matrix_store', nprim,jmax_i_max)
    !allocation 3 
      allocate(matinfo%vec_Rij(3,jmax_i_max,nprim), STAT=istat)
      if(istat .NE. 0) call cq_abort('Allocation 3 in set_matrix_store', nprim,jmax_i_max)
     
   !!set up : idglob_j, beta_j, vec_Rij
    matinfo%idglob_j(:,:) = 0
    matinfo%beta_j(:,:) = 0
    iprim = 0
    ibeg = 1
    matinfo%matrix_size = 0
    do np = 1, bundle%groups_on_node
      if (bundle%nm_nodgroup(np).GT.0) then
        do ni = 1, bundle%nm_nodgroup(np)
          atom_num=bundle%nm_nodbeg(np)+ni-1
          iprim = iprim + 1
          if(atom_num .NE. iprim) write(*,*) 'Warning!! in set_matrix_store : iprim & atom_num = ',iprim,atom_num

          do neigh = 1, mat(np,range)%n_nab(ni)
            ist = mat(np,range)%i_acc(ni) + neigh - 1  ! Neighbour-labeling
            gcspart=BCS_parts%icover_ibeg(mat(np,range)%i_part(ist))+mat(np,range)%i_seq(ist)-1
            j_global_part = BCS_parts%lab_cell(mat(np,range)%i_part(ist))

            matinfo%idglob_j(neigh,iprim) = id_glob(parts%icell_beg(j_global_part) &
                                                +mat(np,range)%i_seq(ist)-1)
            matinfo%beta_j(neigh,iprim) = nsf_species(species_glob(matinfo%idglob_j(neigh,iprim)))
            matinfo%vec_Rij(1,neigh,iprim) = bundle%xprim(iprim) - BCS_parts%xcover(gcspart)
            matinfo%vec_Rij(2,neigh,iprim) = bundle%yprim(iprim) - BCS_parts%ycover(gcspart)
            matinfo%vec_Rij(3,neigh,iprim) = bundle%zprim(iprim) - BCS_parts%zcover(gcspart)

          enddo !neigh = 1, mat(np,range)%n_nab(ni)

          len = matinfo%jbeta_max_i(iprim)*nsf_species(bundle%species(iprim))
          matinfo%matrix_size = matinfo%matrix_size + len
          matinfo%ibeg_data_matrix(iprim) = ibeg
          ibeg = ibeg + len
        enddo !(ni)
      endif !(bundle%nm_nodgroup(np).GT.0) then
    enddo !(np)

    allocate(matinfo%data_matrix(matinfo%matrix_size), STAT = istat)
     if(istat .NE. 0) call cq_abort('Allocation 4 in set_matrix_store', istat, matinfo%matrix_size)
    matinfo%data_matrix(1:matinfo%matrix_size)=mat_p(matA)%matrix(1:matinfo%matrix_size)

    return
   end subroutine set_matrix_store

   subroutine free_matrix_store(matinfo)
    use GenComms, ONLY: cq_abort
    implicit none
    type(matrix_store), intent(inout) :: matinfo
    integer :: istat
   
     deallocate(matinfo%data_matrix, matinfo%vec_Rij, matinfo%beta_j, matinfo%idglob_j, &
                matinfo%ibeg_data_matrix, matinfo%jbeta_max_i, matinfo%jmax_i, matinfo%idglob_i, &
                matinfo%nsf_spec_i, STAT=istat)
     if(istat .NE. 0) call cq_abort("Error in deallocation at free_matrix_store",istat)
    return
   end subroutine free_matrix_store

end module store_matrix
